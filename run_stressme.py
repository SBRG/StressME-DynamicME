import yaml
import pickle
import numpy as np
import sys
from pathlib import Path
import shutil

import cobrame
import ecolime
import AcidifyME
from AcidifyME import *
from AcidifyME.periplasmic_proteome import *
from AcidifyME.membrane_lipid_constraint import *
from AcidifyME.membrane_protein_activity import *
from AcidifyME.proton_influx import *
from oxidizeme.model import *
from qminospy.me1 import ME_NLP1

def main(config_path):
    # --- Load config ---
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    # --- Load ME model ---
    with open(cfg["model_file"], "rb") as f:
        me = pickle.load(f)

    # --- Optionally update keff for DXPRIi and others ---
    for rxn_id, keff_val in cfg.get("updated_dxr_keff", {}).items():
        try:
            rxn = me.reactions.get_by_id(rxn_id)
            rxn.keff = keff_val
            rxn.update()
            print(f"Updated keff for {rxn_id} to {keff_val}")
        except Exception as e:
            print(f"Could not update keff for {rxn_id}: {e}")

    # --- Set parameters ---
    temperature = cfg.get("temperature", 37)
    ROS = cfg.get("ROS", 1.0)
    pH = cfg.get("pH", 7.0)
    me.unmodeled_protein_fraction = cfg.get("unmodeled_protein_fraction", 0.1)

    # --- Always apply thermal stress ---
    ecolime.chaperones.change_temperature(me, temperature)

    # --- Always apply acid stress ---
    fold_change = (
        5.666666666669649 * pH ** 4
        - 1.396666666667337e+02 * pH ** 3
        + 1.286583333333885e+03 * pH ** 2
        - 5.254083333335287e+03 * pH
        + 8.037000000002481e+03
    )
    print(f"fold_change: {fold_change}")
    temperatureK = 273.15 + temperature
    add_membrane_constraint(me, 'MJR_habituated')
    add_protein_stability_constraint(me, pH, temperatureK, constrain_periplasm_biomass=False)
    add_proton_leak_rxn(me, pH)
    modify_membrane_protein_activity(me, fold_change)


    solver = ME_NLP1(me)

    # --- Apply oxidative stress ---
    if ROS > 0:
        stress = StressME(solver)
        stress.make_stressme()
        subs_dict = {}
        for k, v in cfg.get("ros_substrate_dict", {'h2o2_c': 5e-08, 'o2s_c': 2e-10}).items():
            subs_dict[k] = float(v) * ROS
        stress.substitute_ros(solver, subs_dict)
        stress.substitute_metal(solver)

        # Unlock ROS detox enzymes
        for r in me.process_data.SPODM.parent_reactions:
            r.upper_bound = 1000
        for r in me.process_data.SPODMpp.parent_reactions:
            r.upper_bound = 1000
        for r in me.process_data.CAT.parent_reactions:
            r.upper_bound = 1000
        me.reactions.MOX_REV_CPLX_dummy.upper_bound = 0

        stress.force_o2s_from_pq()
        for r in me.reactions.query('PQ2RED'):
            print(r.id, '\t', r.lower_bound)

    # --- Set exchange reaction bounds from config ---
    for rxn_id, lb in cfg.get("exchange_bounds", {}).items():
        try:
            rxn = me.reactions.get_by_id(rxn_id)
            rxn.lower_bound = float(lb)
            print(f"Set lower_bound for {rxn_id} to {lb}")
        except Exception as e:
            print(f"Could not set lower_bound for {rxn_id}: {e}")

    # --- Run simulation (warm start if basis.npy exists) ---
    if cfg.get("warm_start", True):
        try:
            saved_basis = np.load(cfg.get("basis_file", "basis.npy"))
            MU_PREC = float(cfg.get("mu_precision", 1e-4))
            sol, hs, xopt, cache = solver.bisectmu(precision=MU_PREC, basis=saved_basis)
        except Exception as e:
            print("Warm start failed:", e)
            sol, hs, xopt, cache = solver.bisectmu(precision=MU_PREC)
    else:
        MU_PREC = float(cfg.get("mu_precision", 1e-4))
        sol, hs, xopt, cache = solver.bisectmu(precision=MU_PREC)

    # --- Results ---
    print("Solution status:", me.solution.status)
    print("Growth rate:", me.reactions.biomass_dilution.x)

    # --- Save solution ---
    results_dir = Path(cfg['project_root']) / cfg['project_results_folder'] / cfg['project_name']
    results_dir.mkdir(parents=True, exist_ok=True)

    # Save solution inside results_dir
    filename = results_dir / (
        cfg.get("output_prefix", "StressME_heatevolved") +
        f"_T_{temperature}_pH_{pH}_ROS_{ROS}X_sol.pickle"
    )
    with open(filename, "wb") as f:
        pickle.dump(me.solution, f)
    print(f"Solution saved to {filename}")

    # Copy config file to results_dir
    config_yaml_path = results_dir / f"{cfg['project_name']}_config.yaml"
    shutil.copy(config_path, config_yaml_path)
    print(f"→ Config copied to {config_yaml_path}")

    # Save updated ME model inside results_dir
    updated_model_path = results_dir / f"{cfg['project_name']}_updated_model.pickle"
    with open(updated_model_path, "wb") as f:
        pickle.dump(me, f)
    print(f"→ Updated ME model saved to {updated_model_path}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 run_stressme.py <config_stressME.yaml>")
        sys.exit(1)
    main(sys.argv[1])