import pandas as pd
import json
import matplotlib.pyplot as plt
import yaml, pickle
from pathlib import Path
import numpy as np
import uuid
from collections import defaultdict
from IPython.display import display, HTML

cog_colors = {
#     same as https://help.ezbiocloud.net/cog-colors/
    "J": {"color" : "#ff0000", "desc":   "Translation, ribosomal structure and biogenesis"}, # (red)
    "A": {"color" : "#c2af58",  "desc":  "RNA processing and modification"}, # (blue)
    "K": {"color" : "#ff9900",  "desc":  "Transcription"}, # (green)
    "L": {"color" : "#ffff00",  "desc":  "Replication, recombination and repair"}, # (purple)
    "B": {"color" : "#ffc600",  "desc":  "Chromatin structure and dynamics"}, # (orange)
    "D": {"color" : "#99ff00",  "desc":  "Cell cycle control/division"}, # (yellow)
    "Y": {"color" : "#493126",  "desc":  "Nuclear structure"}, # (brown)
    "V": {"color" : "#ff008a",  "desc":  "Defense mechanisms"}, # (pink)
    "T": {"color" : "#0000ff",  "desc":  "Signal transduction mechanisms"}, # (grey)
    "M": {"color" : "#9ec928",  "desc":  "Cell wall/membrane/envelope biogenesis"}, # (teal)
    "N": {"color" : "#006633",  "desc":  "Cell motility"}, # (salmon)
    "Z": {"color" : "#660099",  "desc":  "Cytoskeleton"}, # (lavender-blue)
    "W": {"color" : "#336699",  "desc":  "Extracellular structures"}, # (pink-purple)
    "U": {"color" : "#33cc99",  "desc":  "Intracellular trafficking"}, # (lime green)
    "O": {"color" : "#00ffff",  "desc":  "Posttranslational modification/chaperones"}, # (gold)
    "C": {"color" : "#9900ff",  "desc":  "Energy production and conversion"}, # (tan)
    "G": {"color" : "#805642",  "desc":  "Carbohydrate metabolism"}, # (light grey)
    "E": {"color" : "#ff00ff",  "desc":  "Amino acid metabolism"}, # (dark green)
    "F": {"color" : "#99334d",  "desc":  "Nucleotide metabolism"}, # (orange-red)
    "H": {"color" : "#727dcc",  "desc":  "Coenzyme metabolism"}, # (blue-purple)
    "I": {"color" : "#5c5a1b",  "desc":  "Lipid metabolism"}, # (magenta)
    "P": {"color" : "#0099ff",  "desc":  "Inorganic ion metabolism"}, # (green)
    "Q": {"color" : "#ffcc99",  "desc":  "Secondary metabolites"}, # (mustard)
    "R": {"color" : "#ff9999",  "desc":  "General function prediction"}, # (brown)
    "S": {"color" : "#808080",  "desc":  "Function unknown"}, # (dark grey),
    "Unknown": { "color" : "#bbbbbb", "desc" : "Unknown"} #(gray)
}

# --- loaders ---
def load_config_and_model(config_path):
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)
    with open(cfg["model_file"], "rb") as f:
        me = pickle.load(f)
    return cfg, me

def load_results(results_path):
    # Load results from JSON
    with open(results_path, 'r') as f:
        results = json.load(f)
    return results

def mw_dict_from_me(me):
    # g / mmol (kDa == g/mmol). Some ME objects store 'mass'; others 'formula_weight'
    return {m.id: getattr(m, "mass", getattr(m, "formula_weight", np.nan)) for m in me.metabolites}

# --- unit conversion ---
def concentrations_to_df(concentrations):
    # concentrations: list[dict] from DynamicME result["concentration"]
    return pd.DataFrame(concentrations)

def convert_mM_to_gL(df_mM, mw_dict):
    df = df_mM.copy()
    for col in df.columns:
        if col in mw_dict and pd.notna(mw_dict[col]):
            df[col] = df[col] * mw_dict[col] / 1000.0  # (mM)*(g/mmol)/1000 => g/L
    return df

# --- plotting ---
def plot_concentrations(time, df_conc, yunit="mM"):
    df = df_conc.copy()
    df["time"] = time
    plt.figure(figsize=(8,5))
    for met in df.columns:
        if met != "time":
            plt.plot(df["time"], df[met], marker="o", label=met)
    plt.xlabel("Time (hours)")
    plt.ylabel(f"Concentration ({yunit})")
    plt.title("Metabolite Concentrations Over Time")
    plt.legend(loc=(1.05, 0), title="Tracked Metabolites")
    plt.grid(True)
    plt.show()

def plot_growth_and_mu(time, biomass, rxn_flux_list, growth_rxn_id="biomass_dilution"):
    mu = np.array([rf.get(growth_rxn_id, np.nan) for rf in rxn_flux_list])
    fig, ax1 = plt.subplots(figsize=(6,4))

    # Biomass in blue
    color1 = "tab:blue"
    ax1.plot(time, biomass, marker="o", color=color1, label="Biomass")
    ax1.set_xlabel("Time (hours)")
    ax1.set_ylabel("Biomass (gDW/L)", color=color1)
    ax1.tick_params(axis="y", labelcolor=color1)

    # Growth rate μ in red
    ax2 = ax1.twinx()
    color2 = "tab:red"
    ax2.plot(time, mu, marker="s", linestyle="--", color=color2, label="Growth rate (μ)")
    ax2.set_ylabel("Growth rate μ (hr⁻¹)", color=color2)
    ax2.tick_params(axis="y", labelcolor=color2)

    ax1.grid(True)

    # Combine legends
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, loc=(1.05,0))

    plt.title("Biomass & Growth Rate")
    plt.show()
    
def plot_translation_fluxes(results, cog_df=None, cog_colors=None, top_k=None):
    """
    Plots translation reaction trajectories two ways:
      1) Flux normalized by biomass (mmol gDW^-1 hr^-1)
      2) Absolute synthesis rate (mmol L^-1 hr^-1) = flux * biomass

    Args:
        results: dict from DynamicME simulate_batch
        cog_df: DataFrame with COG annotations, indexed by locus (or gene).
        cog_colors: dict mapping COG category letters -> color.
        top_k: if set, plot only the top_k reactions by max absolute rate.
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    time = np.asarray(results["time"])
    biomass = np.asarray(results["biomass"])  # gDW/L
    df_flux = pd.DataFrame(results["rxn_flux"])   # per-gDW fluxes

    # pick translation reactions that are ever > 0
    tx_cols = [c for c in df_flux.columns
               if c.startswith("translation_") and df_flux[c].max() > 0]

    if not tx_cols:
        print("No translation reaction fluxes > 0 found.")
        return

    # per-gDW flux
    df_flux_tx = df_flux[tx_cols].fillna(0.0)
    # absolute rate (flux * biomass)
    df_rate_abs = df_flux_tx.mul(biomass, axis=0).fillna(0.0)

    # optional: restrict to top_k
    if top_k is not None and top_k < len(tx_cols):
        peak = df_rate_abs.max(axis=0).sort_values(ascending=False)
        tx_cols = list(peak.index[:top_k])
        df_flux_tx = df_flux_tx[tx_cols]
        df_rate_abs = df_rate_abs[tx_cols]

    # --- plotting ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    for col in tx_cols:
        gene = col.replace("translation_", "")

        # default color
        color = None
        if cog_df is not None and cog_colors is not None:
            # lookup by gene or locus
            if gene in cog_df.index:
                cat = str(cog_df.loc[gene, "COG category"])
                color = get_cog_color(cat, cog_colors)
            elif "gene" in cog_df.columns and gene in cog_df["gene"].values:
                row = cog_df[cog_df["gene"] == gene].iloc[0]
                cat = str(row["COG category"])
                color = get_cog_color(cat, cog_colors)

        ax1.plot(time, df_flux_tx[col], label=gene, color=color)
        ax2.plot(time, df_rate_abs[col], label=gene, color=color)

    ax1.set_ylabel("Translation Flux (mmol gDW$^{-1}$ hr$^{-1}$)")
    ax1.set_title("Translation Reactions – per biomass")
    ax1.grid(True)

    ax2.set_xlabel("Time (hours)")
    ax2.set_ylabel("Synthesis Rate (mmol L$^{-1}$ hr$^{-1}$)")
    ax2.set_title("Translation Reactions – absolute (per liter)")
    ax2.grid(True)

    handles, labels = ax2.get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()

def plot_formation_fluxes(results, top_k=None, show_absolute=True):
    """
    Plot protein complex 'formation_' reaction trajectories.

    Panels:
      (A) Flux normalized by biomass (mmol gDW^-1 hr^-1)
      (B) Absolute synthesis rate (mmol L^-1 hr^-1) = flux * biomass  [if show_absolute]

    Args:
        results: dict returned by DynamicME simulate_batch
        top_k:   if set, only the top_k reactions by peak absolute rate are shown
        show_absolute: include the absolute-rate panel
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    time = np.asarray(results["time"])
    biomass = np.asarray(results["biomass"])  # gDW/L
    df_flux = pd.DataFrame(results["rxn_flux"]).fillna(0.0)

    # pick formation_* reactions that are ever > 0
    form_cols = [c for c in df_flux.columns if c.startswith("formation_") and df_flux[c].max() > 0]
    if not form_cols:
        print("No formation reaction fluxes > 0 found.")
        return

    df_form = df_flux[form_cols]

    # absolute rate (mmol L^-1 hr^-1) = flux * biomass
    df_abs = df_form.mul(biomass, axis=0)

    # optional: keep only strongest top_k by peak absolute rate
    if top_k is not None and top_k < len(form_cols):
        peak = df_abs.max(axis=0).sort_values(ascending=False)
        form_cols = list(peak.index[:top_k])
        df_form = df_form[form_cols]
        df_abs  = df_abs[form_cols]

    # set up figure
    if show_absolute:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    else:
        fig, ax1 = plt.subplots(1, 1, figsize=(8, 4))
        ax2 = None

    # (A) per‑gDW flux
    for col in form_cols:
        ax1.plot(time, df_form[col], label=col)
    ax1.set_ylabel("Formation Flux (mmol gDW$^{-1}$ hr$^{-1}$)")
    ax1.set_title("Formation Reaction Fluxes – per biomass")
    ax1.grid(True)
    ax1.set_ylim(bottom=0)

    # (B) absolute per‑liter rate (optional)
    if show_absolute:
        for col in form_cols:
            ax2.plot(time, df_abs[col], label=col)
        ax2.set_xlabel("Time (hours)")
        ax2.set_ylabel("Formation Rate (mmol L$^{-1}$ hr$^{-1}$)")
        ax2.set_title("Formation Reaction Rates – absolute (per liter)")
        ax2.grid(True)
        ax2.set_ylim(bottom=0)
        handles, labels = ax2.get_legend_handles_labels()
    else:
        ax1.set_xlabel("Time (hours)")
        handles, labels = ax1.get_legend_handles_labels()

    # optional outside legend (uncomment to show)
    # fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1.02, 0.5))

    fig.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()
    
def plot_biomass_dilution_fluxes(results, top_k=None, show_absolute=True):
    """
    Plot biomass_to_biomass reaction fluxes with clean legend labels.
    """
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    time = np.asarray(results["time"])
    biomass = np.asarray(results["biomass"])  # gDW/L
    df_flux = pd.DataFrame(results["rxn_flux"]).fillna(0.0)

    # pick biomass_to_biomass reactions
    rxn_cols = [c for c in df_flux.columns if "biomass_to_biomass" in c and df_flux[c].max() > 0]
    if not rxn_cols:
        print("No biomass_to_biomass reaction fluxes > 0 found.")
        return

    df_sel = df_flux[rxn_cols]
    df_abs = df_sel.mul(biomass, axis=0)  # mmol L^-1 hr^-1

    # optionally keep only strongest reactions
    if top_k is not None and top_k < len(rxn_cols):
        peak = df_abs.max(axis=0).sort_values(ascending=False)
        rxn_cols = list(peak.index[:top_k])
        df_sel = df_sel[rxn_cols]
        df_abs = df_abs[rxn_cols]

    # make pretty labels: strip "_biomass_to_biomass"
    def pretty_label(rxn_id):
        return rxn_id.replace("_biomass_to_biomass", "").replace("_", " ").title()

    labels = [pretty_label(c) for c in rxn_cols]

    if show_absolute:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    else:
        fig, ax1 = plt.subplots(1, 1, figsize=(8, 4))
        ax2 = None

    # per-gDW flux
    for col, lbl in zip(rxn_cols, labels):
        ax1.plot(time, df_sel[col], label=lbl)
    ax1.set_ylabel("Flux (mmol gDW$^{-1}$ hr$^{-1}$)")
    ax1.set_title("Biomass Component Fluxes – per biomass")
    ax1.grid(True)

    # absolute per-liter rate
    if show_absolute:
        for col, lbl in zip(rxn_cols, labels):
            ax2.plot(time, df_abs[col], label=lbl)
        ax2.set_xlabel("Time (hours)")
        ax2.set_ylabel("Rate (mmol L$^{-1}$ hr$^{-1}$)")
        ax2.set_title("Biomass Component Rates – absolute (per liter)")
        ax2.grid(True)
        handles, labels = ax2.get_legend_handles_labels()
    else:
        ax1.set_xlabel("Time (hours)")
        handles, labels = ax1.get_legend_handles_labels()

    # put legend outside
    fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1.15, 0.5), title="Biomass Components")
    fig.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()



def export_translation_weights(me,results, cog_df=None, t_query=None,
                               value="flux", out_file=None):
    """
    Export translation reaction weights at a given time point.

    Args:
        results: dict from DynamicME simulate_batch
        cog_df: DataFrame with COG annotations (optional, for locus/gene mapping)
        t_query: float, time in hours (closest snapshot). If None, take last time point.
        value: "flux" (mmol gDW^-1 hr^-1), "rate" (mmol L^-1 hr^-1), or "mass" (g L^-1 hr^-1 approx).
        out_file: optional path to write TSV for ProteoMaps ("locus<TAB>weight")

    Returns:
        DataFrame with columns [locus, gene, COG category, weight]
    """
    time = np.asarray(results["time"])
    biomass = np.asarray(results["biomass"])
    df_flux = pd.DataFrame(results["rxn_flux"])

    # pick translation reactions
    tx_cols = [c for c in df_flux.columns if c.startswith("translation_")]
    if not tx_cols:
        raise ValueError("No translation reactions found in results.")

    # pick index
    if t_query is None:
        ti = -1
    else:
        ti = int(np.argmin(np.abs(time - t_query)))

    flux_row = df_flux.loc[ti, tx_cols].fillna(0.0)  # mmol gDW^-1 hr^-1
    biomass_t = biomass[ti]

    rows = []
    for rxn_id, flux in flux_row.items():
        if flux <= 0:
            continue
        gene = rxn_id.replace("translation_", "")
        p = me.metabolites.get_by_id('protein_{}'.format(gene))
        locus = gene
        cat = None
        if cog_df is not None:
            try:
                # safest check: ensure exactly one row is matched
                if gene in cog_df.index:
                    row = cog_df.loc[[gene]].iloc[0]  # get the first row even if duplicate index
                    locus = gene
                    cat = str(row.get("COG category", ""))
                    gene_name = row.get("gene", gene)
                elif "gene" in cog_df.columns and gene in cog_df["gene"].values:
                    row = cog_df[cog_df["gene"] == gene].iloc[0]
                    locus = row.name if row.name else gene
                    cat = str(row.get("COG category", ""))
                    gene_name = row.get("gene", gene)
                else:
                    gene_name = gene
            except Exception as e:
                print(f"⚠️ Failed to map gene '{gene}': {e}")
                gene_name = gene
                cat = ""


        if value == "flux":
            weight = flux  # mmol gDW^-1 hr^-1
        elif value == "rate":
            weight = flux * biomass_t  # mmol L^-1 hr^-1
        elif value == "mass":
            # approximate: rate * MW (g/mmol)
            weight = flux * biomass_t * (p.formula_weight / 1000.0)

        else:
            raise ValueError("value must be 'flux', 'rate', or 'mass'.")

        rows.append({
            "locus": locus,
            "gene": gene_name,
            "COG category": cat,
            "weight": weight
        })

    df_out = pd.DataFrame(rows)
    if out_file:
        df_out[["locus", "weight"]].to_csv(out_file, sep="\t", index=False)
        print(f"Wrote {len(df_out)} rows to {out_file}")

    return df_out.sort_values("weight", ascending=False)



# === Step 1: Color assignment function ===
# def get_cog_color(cat, cog_colors):
#     if not isinstance(cat, str):
#         return cog_colors.get("Unknown", "#bbbbbb")

#     cat = cat.strip()
#     if cat == "" or cat.lower() in {"none", "unknown"}:
#         return cog_colors.get("Unknown", "#bbbbbb")

#     for letter in cat:
#         if letter in cog_colors:
#             return cog_colors[letter]

#     return cog_colors.get("Unknown", "#bbbbbb")

def get_cog_color(cat, cog_colors):
    if not isinstance(cat, str):
        return cog_colors.get("Unknown", {}).get("color", "#bbbbbb")

    cat = cat.strip()
    if cat == "" or cat.lower() in {"none", "unknown"}:
        return cog_colors.get("Unknown", {}).get("color", "#bbbbbb")

    for letter in cat:
        if letter in cog_colors:
            return cog_colors[letter]["color"]

    return cog_colors.get("Unknown", {}).get("color", "#bbbbbb")



# === Step 2: Build hierarchy ===
def build_hierarchy_for_d3_from_df(df, cog_col="COG category", name_col="locus", weight_col="weight", cog_colors=None, label_col="locus"):
    cat_to_genes = defaultdict(list)

    for _, row in df.iterrows():
        cat = row[cog_col] if pd.notna(row[cog_col]) else "Unknown"
        name = row[label_col]
        weight = float(row[weight_col])
        color = get_cog_color(cat, cog_colors)
        if weight > 0:
            cat_to_genes[cat].append({"name": name, "weight": weight, "color": color})

    children = [{"name": cat, "children": genes} for cat, genes in cat_to_genes.items()]
    return {"name": "root", "children": children}


# === Step 3: Write interactive HTML Voronoi treemap ===
def write_d3_voronoi_html(hierarchy, out_html="voronoi.html", width=800, height=800, seed="qsproteome", cog_colors=None):
    import json, uuid
    from IPython.display import display, HTML

    uid = uuid.uuid4().hex

    # Sanitize cog_colors into JS-compatible string
    cog_colors_js = json.dumps(cog_colors or {}, indent=2)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Voronoi Treemap</title>
  <style>
    .cell {{ stroke: #fff; stroke-width: 0.5; }}
    text {{
      pointer-events: none;
      text-anchor: middle;
    }}
  </style>
  <script src="https://d3js.org/d3.v6.min.js"></script>
  <script src="https://rawcdn.githack.com/Kcnarf/d3-weighted-voronoi/v1.1.3/build/d3-weighted-voronoi.js"></script>
  <script src="https://rawcdn.githack.com/Kcnarf/d3-voronoi-map/v2.1.1/build/d3-voronoi-map.js"></script>
  <script src="https://rawcdn.githack.com/Kcnarf/d3-voronoi-treemap/v1.1.2/build/d3-voronoi-treemap.js"></script>
</head>
<body>
<div style="display: flex;">
  <svg id="svg-{uid}" width="{width}" height="{height}"></svg>
  <div style="margin-left: 20px;">
    <h3>COG Categories</h3>
    <ul style="list-style: none; padding: 0;" id="legend-{uid}"></ul>
  </div>
</div>

<script>
  const width = {width}, height = {height};
  const data = {json.dumps(hierarchy)};
  const root = d3.hierarchy(data).sum(d => d.weight || 0);
  const treemap = d3.voronoiTreemap()
    .clip([[0, 0], [0, height], [width, height], [width, 0]])
    .prng(d3.randomLcg("{seed}"));

  treemap(root);

  const svg = d3.select("#svg-{uid}");

  // Draw cells
  svg.selectAll("path")
    .data(root.leaves())
    .enter()
    .append("path")
    .attr("d", d => d3.line()(d.polygon) + "Z")
    .attr("class", "cell")
    .attr("fill", d => d.data.color || "#cccccc");

  // Compute font scale
  const weights = root.leaves().map(d => d.data.weight || 0);
  const minWeight = d3.min(weights);
  const maxWeight = d3.max(weights);
  const fontSizeScale = d3.scaleSqrt()
    .domain([minWeight, maxWeight])
    .range([10, 30]);

  // Add labels
  svg.selectAll("text")
    .data(root.leaves())
    .enter()
    .append("text")
    .attr("x", d => d3.polygonCentroid(d.polygon)[0])
    .attr("y", d => d3.polygonCentroid(d.polygon)[1])
    .attr("font-size", d => fontSizeScale(d.data.weight))
    .text(d => d.data.name);

  // Legend
  const cogColors = {cog_colors_js};
  const legendContainer = d3.select("#legend-{uid}");
  Object.entries(cogColors).forEach(([letter, info]) => {{
    const item = legendContainer.append("li").style("margin-bottom", "4px");
    item.append("span")
      .style("display", "inline-block")
      .style("width", "12px")
      .style("height", "12px")
      .style("background-color", info.color)
      .style("margin-right", "6px");
    item.append("span").text(`${{letter}}: ${{info.desc}}`);
  }});
</script>
</body>
</html>"""

    with open(out_html, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"✅ Wrote: {out_html}")
    display(HTML(f'<iframe src="{out_html}" width="{width+300}" height="{height}"></iframe>'))
