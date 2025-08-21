import yaml
import sys
from pathlib import Path
import os
import pickle
import pandas as pd
import numpy as np
import json
import shutil


def main(config_path):
    # Load config
    print ("→ Loading config file from...", config_path)

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Prepare output paths
    results_dir = Path(config['project_root']) / config['project_results_folder'] / config['project_name']
    results_dir.mkdir(parents=True, exist_ok=True)
    results_json_path = results_dir / f"{config['project_name']}_results.json"
    config_yaml_path = results_dir / f"{config['project_name']}_config.yaml"
    print(f"→ Results directory: {results_dir}")
    print(f"→ Results JSON path: {results_json_path}")
    print(f"→ Config YAML path: {config_yaml_path}")

    # Set up paths
    project_root = Path(config.get('project_root', Path.cwd()))
    sys.path.insert(0, str(project_root))

    # Load model
    print ("→ Loading model file from...", config['model_file'])
    with open(config['model_file'], 'rb') as f:
        me = pickle.load(f)

    # Import DynamicME and ParamOpt
    from dynamicme.dynamic import DynamicME, ParamOpt
    from cobrame.core.reaction import MetabolicReaction

    print ("→ Loading simulation parameters:")
    print(f"   dt    : {config['dt']}")
    print(f"   T     : {config['T']}")
    print(f"   V     : {config['V']}")
    print(f"   X0    : {config['X0']}")
    print(f"   LB_EX : {config['LB_EX']}")
    print(f"   LB_O2 : {config['LB_O2']}")
    # Simulation parameters
    T = config['T']
    dt = config['dt']
    V = config['V']
    X0 = config['X0']
    c0_dict = config['c0_dict']
    LB_EX = config['LB_EX']
    LB_O2 = config['LB_O2']


    extra_rxns_tracked = [me.reactions.biomass_dilution, me.reactions.EX_o2_e]

    # Load tracked metabolites from config
    tracked_met_ids = config.get('tracked_metabolites', [])
    tracked_mets = [me.metabolites.get_by_id(mid) for mid in tracked_met_ids]
    print("→ Tracking metabolites:")
    for met in tracked_mets:
        print(f"   {met.id:10} : {met.name}")

    # Track metabolic reactions for each metabolite
    print("→ Tracking metabolic reactions:")
    rows_tracked = []
    for met_track in tracked_mets:
        mid_c = met_track.id.replace('_p','_c')
        mid_e = met_track.id.replace('_p','_e')
        met_c = me.metabolites.get_by_id(mid_c)
        met_e = me.metabolites.get_by_id(mid_e)
        for rxn in met_track.reactions:
            if isinstance(rxn, MetabolicReaction) and rxn.keff and (met_c in rxn.metabolites or met_e in rxn.metabolites):
                extra_rxns_tracked.append(rxn)
                rows_tracked.append({'met':met_track.id, 'rxn':rxn.id})
                print(f"   {rxn.id:25} : {rxn.reaction:50}")


    df_tracked = pd.DataFrame(rows_tracked)

    # Track translation fluxes
    tracked_trsl_cfg = config.get('tracked_translation_reactions', False)
    if tracked_trsl_cfg is True:
        rxns_trsl = me.reactions.query('translation_')
        print(f"→ Tracking {len(rxns_trsl)} translations reactions from ME model...")
        print(f"   example:")
        print(f"   {rxns_trsl[0].id:25} : {str(rxn.reaction)[0:200]}...")
    elif isinstance(tracked_trsl_cfg, list):
        rxns_trsl = [me.reactions.get_by_id(rid) for rid in tracked_trsl_cfg]
        print(f"→ Tracking {len(rxns_trsl)} translations reactions from config...")
        for rxn in rxns_trsl:
            print(f"   {rxn.id:25} : {str(rxn.reaction)[0:200]}...")
    else:
        rxns_trsl = []

    extra_rxns_tracked += rxns_trsl
                

    # Track biomass reactions
    tracked_biomass_cfg = config.get('tracked_biomass_to_biomass_reactions', False)
    if tracked_biomass_cfg is True:
        biomass_rxns = [r for r in me.reactions if r.id.endswith('_biomass_to_biomass')]
        print(f"→ Tracking {len(biomass_rxns)} biomass reactions from ME model...")
    elif isinstance(tracked_biomass_cfg, list):
        biomass_rxns = [me.reactions.get_by_id(rid) for rid in tracked_biomass_cfg]
        print(f"→ Tracking  {len(biomass_rxns)} biomass reactions from config...")
    else:
        biomass_rxns = []
    extra_rxns_tracked += biomass_rxns
    for rxn in biomass_rxns:
        print(f"   {rxn.id:25} : {rxn.reaction:50}")

     # Track complex formation reactions
    tracked_complex_cfg = config.get('tracked_complex_formation_reactions', False)
    print ("tracked_complex_cfg",tracked_complex_cfg)
    if tracked_complex_cfg is True:
        complex_rxns = [r for r in me.reactions if r.id.startswith('formation_')]
        print(f"→ Tracking {len(complex_rxns)} complex formation reactions from ME model...")
    elif isinstance(tracked_complex_cfg, list):
        complex_rxns = [me.reactions.get_by_id(rid) for rid in tracked_complex_cfg]
        print(f"→ Tracking  {len(complex_rxns)} complex formation reactions from config...")
    else:
        complex_rxns = []
    extra_rxns_tracked += complex_rxns
    
    # Convert concentrations
    print ("→ Setting up initial concentrations from config...")
    for mid, c in c0_dict.items():
        met = me.metabolites.get_by_id(mid)
        c0_dict[met.id] = c / met.formula_weight * 1000

    # Set max uptake rates
    print ("→ Setting up max uptake rates...")
    lb_dict = {}
    for mid in c0_dict.keys():
        rxn = DynamicME(me).get_exchange_rxn(mid)
        lb_dict[rxn.id] = LB_O2 if rxn.id == 'EX_o2_e' else LB_EX
    lb_dict[me.reactions.EX_o2_e.id] = LB_O2
    ub_dict = {}

    # Run simulation
    dyme = DynamicME(me)
    # return
    print ("Starting DynamicME simulation...")

    result = dyme.simulate_batch(
        T=T, #hours
        dt = dt, #time-step (hours)
        c0_dict=c0_dict,
        X0=X0 / V, #g/L
        prec_bs=1e-3,
        ZERO_CONC=0.,
        extra_rxns_tracked=extra_rxns_tracked,
        lb_dict=lb_dict,
        ub_dict=ub_dict,
        verbosity=1,
    )


    # Save results as JSON
    def convert_for_json(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, dict):
            return {k: convert_for_json(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [convert_for_json(v) for v in obj]
        return obj

    json_result = convert_for_json(result)
    with open(results_json_path, 'w') as f:
        json.dump(json_result, f)
    print(f"→ Results saved to {results_json_path}")

    # Save config as YAML
    shutil.copy(config_path, config_yaml_path)
    print(f"→ Config copied to {config_yaml_path}")



    # # Convert solution's metabolite concentrations from mM to g/L
    # sim_params = {
    #     'T': T,
    #     'X0': X0,
    #     'c0_dict': c0_dict,
    #     'lb_dict': lb_dict,
    #     'ub_dict': ub_dict,
    #     'extra_rxns_tracked': extra_rxns_tracked,
    #     'ZERO_CONC': 0.
    # }
    # growth_rxn = me.reactions.biomass_dilution
    # popt = ParamOpt(me, sim_params, growth_rxn=growth_rxn.id)
    # sol = popt.compute_conc_profile(result)

    # df_mw = pd.DataFrame([{'id':m.id,'mass':getattr(m,'mass',m.formula_weight),'name':m.name} for m in me.metabolites])
    # sol_gL = sol.copy()
    # variables = c0_dict.keys()
    # for col in sol.columns:
    #     if col in variables:
    #         c_mM = sol[col]
    #         mw = df_mw[df_mw.id==col].mass.values[0]
    #         c_gL = c_mM * mw / 1000.
    #         sol_gL[col] = c_gL

    # # Save results
    # sol_gL.to_csv(config.get('output_file', 'dynamicme_example_output.csv'))
    # pd.DataFrame(result['ex_flux']).to_csv(config.get('ex_flux_file', 'ex_flux.csv'))

    # print(f"Results saved to {config.get('output_file', 'dynamicme_example_output.csv')}")
    # print(f"Exchange fluxes saved to {config.get('ex_flux_file', 'ex_flux.csv')}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 run_dynamicme.py <config_dynamicME.yaml>")
        sys.exit(1)
    main(sys.argv[1])