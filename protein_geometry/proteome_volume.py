"""
Proteome Volume Calculation Module

This module provides functions to calculate the total volume of the proteome
by multiplying translation fluxes from ME model solutions with individual
protein volumes from AlphaFold structures.

Functions:
    load_uniprot_to_blattner_mapping: Parse E. coli gene file for ID mappings
    load_protein_volumes: Load protein volume data from JSONL file
    calculate_proteome_volume_from_solution: Main calculation function
    analyze_proteome_volume_reasonableness: Biological context analysis
    analyze_proteome_volume_by_compartment: Analyze volume distribution by cellular compartment
"""

import json
import re
import numpy as np
from scipy.constants import Avogadro


def load_uniprot_to_blattner_mapping(file_path):
    """
    Parse the E. coli gene file to create mapping from UniProt IDs to Blattner IDs (b-numbers).
    
    Args:
        file_path: Path to the E. coli gene file
        
    Returns:
        dict: {uniprot_id: blattner_id} mapping
    """
    uniprot_to_blattner = {}
    
    with open(file_path, 'r') as f:
        for line_num, line in enumerate(f):
            if line_num == 0:  # Skip header
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                blattner_id = parts[1]  # Accession-1 column (b-number)
                uniprot_column = parts[6]  # UniProt column
                
                # Try multiple patterns to extract UniProt ID
                # Pattern 1: From href attribute (most reliable)
                href_match = re.search(r"uniprot/([PQ][0-9A-Z]+)", uniprot_column)
                if href_match and blattner_id.startswith('b'):
                    uniprot_id = href_match.group(1)
                    uniprot_to_blattner[uniprot_id] = blattner_id
                    continue
                
                # Pattern 2: Direct pattern (fallback)
                direct_match = re.search(r'[PQ][0-9A-Z]+', uniprot_column)
                if direct_match and blattner_id.startswith('b'):
                    uniprot_id = direct_match.group()
                    uniprot_to_blattner[uniprot_id] = blattner_id
    
    return uniprot_to_blattner


def load_protein_volumes(file_path):
    """
    Load protein volumes from all_proteins_results.jsonl
    
    Args:
        file_path: Path to the JSONL file
        
    Returns:
        dict: {uniprot_id: volume_A3} mapping
    """
    volumes = {}
    
    with open(file_path, 'r') as f:
        for line in f:
            data = json.loads(line.strip())
            uniprot_id = data['UniProt_ID']
            volume_A3 = data['vol_A3']  # Volume in Å³
            volumes[uniprot_id] = volume_A3
    
    return volumes


def calculate_proteome_volume_from_solution(me_solution, blattner_to_uniprot, protein_volumes):
    """
    Calculate total proteome volume directly from ME model solution.
    
    The ME solution contains translation fluxes in units of mmol gDW^-1 hr^-1.
    
    Args:
        me_solution: ME model solution object with x_dict attribute
        blattner_to_uniprot: dict mapping Blattner IDs (b-numbers) to UniProt IDs
        protein_volumes: dict mapping UniProt IDs to volumes in Å³
        
    Returns:
        dict: Contains total volume flux and breakdown by protein
    """
    
    solution_dict = me_solution.x_dict
    
    # Find all translation reactions
    translation_fluxes = {k: v for k, v in solution_dict.items() 
                         if k.startswith('translation_') and v > 0}
    
    total_volume_flux_A3_per_gdw_per_hr = 0.0
    protein_breakdown = []
    missing_proteins = []
    
    for rxn_id, flux_mmol_gdw_hr in translation_fluxes.items():
        # Extract Blattner ID from reaction name
        blattner_id = rxn_id.replace("translation_", "")
        
        # Get UniProt ID
        uniprot_id = blattner_to_uniprot.get(blattner_id)
        if not uniprot_id:
            missing_proteins.append(blattner_id)
            continue
            
        # Get protein volume
        volume_A3 = protein_volumes.get(uniprot_id)
        if volume_A3 is None:
            missing_proteins.append(f"{blattner_id} ({uniprot_id})")
            continue
        
        # Convert flux to molecules per gDW per hr
        # flux_mmol_gdw_hr * (1e-3 mol/mmol) * (Avogadro molecules/mol)
        molecules_per_gdw_per_hr = flux_mmol_gdw_hr * 1e-3 * Avogadro
        
        # Calculate volume flux (Å³ per gDW per hr)
        volume_flux_A3_per_gdw_per_hr = molecules_per_gdw_per_hr * volume_A3
        
        total_volume_flux_A3_per_gdw_per_hr += volume_flux_A3_per_gdw_per_hr
        
        protein_breakdown.append({
            'blattner_id': blattner_id,
            'uniprot_id': uniprot_id,
            'flux_mmol_gdw_hr': flux_mmol_gdw_hr,
            'molecules_per_gdw_hr': molecules_per_gdw_per_hr,
            'volume_A3': volume_A3,
            'volume_flux_A3_gdw_hr': volume_flux_A3_per_gdw_per_hr
        })
    
    # Convert to more reasonable units
    # 1 Å³ = 1e-30 m³
    total_volume_flux_m3_per_gdw_per_hr = total_volume_flux_A3_per_gdw_per_hr * 1e-30
    total_volume_flux_uL_per_gdw_per_hr = total_volume_flux_m3_per_gdw_per_hr * 1e9  # to microliters
    
    # Sort by volume contribution
    protein_breakdown.sort(key=lambda x: x['volume_flux_A3_gdw_hr'], reverse=True)
    
    results_dict = {
        'total_volume_flux_A3_per_gdw_per_hr': total_volume_flux_A3_per_gdw_per_hr,
        'total_volume_flux_m3_per_gdw_per_hr': total_volume_flux_m3_per_gdw_per_hr,
        'total_volume_flux_uL_per_gdw_per_hr': total_volume_flux_uL_per_gdw_per_hr,
        'proteins_analyzed': len(protein_breakdown),
        'proteins_missing': len(missing_proteins),
        'missing_protein_list': missing_proteins,
        'protein_breakdown': protein_breakdown,
        'units_explanation': {
            'flux_units': 'mmol gDW^-1 hr^-1 (from ME solution)',
            'volume_units': 'Å³ (from AlphaFold structures)',
            'result_units': 'Å³ gDW^-1 hr^-1 (volume synthesized per gram dry weight per hour)'
        }
    }
    
    return results_dict


def analyze_proteome_volume_reasonableness(proteome_volume, growth_rate_hr):
    """
    Analyze whether the calculated proteome volume flux is biologically reasonable.
    
    Args:
        proteome_volume: dict returned by calculate_proteome_volume_from_solution
        growth_rate_hr: growth rate in hr^-1 from ME solution
        
    Returns:
        dict: Analysis results including assessment and biological context
    """
    
    # E. coli cell properties (typical values)
    ecoli_volume_fL = 0.65  # femtoliters (fL), typical E. coli volume
    ecoli_dry_weight_fg = 280  # femtograms (fg), typical E. coli dry weight
    ecoli_dry_weight_g = ecoli_dry_weight_fg * 1e-15  # convert to grams
    protein_fraction_dry_weight = 0.55  # ~55% of dry weight is protein
    
    # Calculate volume flux per cell
    volume_flux_per_cell_uL_hr = proteome_volume['total_volume_flux_uL_per_gdw_per_hr'] * ecoli_dry_weight_g
    volume_flux_per_cell_fL_hr = volume_flux_per_cell_uL_hr * 1e15
    
    # Compare to cell volume
    volume_flux_as_fraction_of_cell = volume_flux_per_cell_fL_hr / ecoli_volume_fL
    volume_flux_percentage = volume_flux_as_fraction_of_cell * 100
    
    # Expected protein synthesis for growth
    doubling_time_hr = np.log(2) / growth_rate_hr
    protein_volume_doubling_rate_fL_hr = (protein_fraction_dry_weight * ecoli_volume_fL) / doubling_time_hr
    
    # Compare calculated vs expected
    ratio = volume_flux_per_cell_fL_hr / protein_volume_doubling_rate_fL_hr
    
    # Assessment
    if 0.5 <= ratio <= 2.0:
        assessment = "✅ REASONABLE"
        explanation = "The calculated flux is within expected biological range"
    elif 0.1 <= ratio <= 5.0:
        assessment = "⚠️  POSSIBLY REASONABLE"
        explanation = "The calculated flux is somewhat higher/lower than expected but could be valid"
    else:
        assessment = "❌ LIKELY UNREASONABLE"
        explanation = "The calculated flux is significantly different from biological expectations"
    
    return {
        'assessment': assessment,
        'explanation': explanation,
        'ratio_calculated_vs_expected': ratio,
        'volume_flux_per_cell_fL_hr': volume_flux_per_cell_fL_hr,
        'volume_flux_percentage_of_cell': volume_flux_percentage,
        'expected_protein_synthesis_fL_hr': protein_volume_doubling_rate_fL_hr,
        'doubling_time_hr': doubling_time_hr,
        'cell_properties': {
            'volume_fL': ecoli_volume_fL,
            'dry_weight_fg': ecoli_dry_weight_fg,
            'protein_fraction': protein_fraction_dry_weight
        }
    }


def print_proteome_volume_summary(proteome_volume, analysis=None):
    """
    Print a comprehensive summary of proteome volume calculations.
    
    Args:
        proteome_volume: dict returned by calculate_proteome_volume_from_solution
        analysis: dict returned by analyze_proteome_volume_reasonableness (optional)
    """
    
    print("=== PROTEOME VOLUME CALCULATION SUMMARY ===")
    print(f"Total proteins analyzed: {proteome_volume['proteins_analyzed']}")
    print(f"Proteins missing volume data: {proteome_volume['proteins_missing']}")
    print()
    
    print("=== TOTAL VOLUME FLUX ===")
    print(f"Total volume flux: {proteome_volume['total_volume_flux_A3_per_gdw_per_hr']:.2e} Å³ gDW⁻¹ hr⁻¹")
    print(f"Total volume flux: {proteome_volume['total_volume_flux_uL_per_gdw_per_hr']:.6f} μL gDW⁻¹ hr⁻¹")
    print()
    
    print("=== TOP 10 PROTEINS BY VOLUME CONTRIBUTION ===")
    for i, protein in enumerate(proteome_volume['protein_breakdown'][:10]):
        pct = (protein['volume_flux_A3_gdw_hr'] / proteome_volume['total_volume_flux_A3_per_gdw_per_hr']) * 100
        print(f"{i+1:2d}. {protein['blattner_id']:6s} ({protein['uniprot_id']:7s}): "
              f"{protein['volume_flux_A3_gdw_hr']:.2e} Å³/gDW/hr ({pct:.1f}%)")
    
    if analysis:
        print()
        print("=== BIOLOGICAL REASONABLENESS ===")
        print(f"{analysis['assessment']}")
        print(f"{analysis['explanation']}")
        print(f"Calculated/Expected ratio: {analysis['ratio_calculated_vs_expected']:.2f}")
        print(f"Volume flux per cell: {analysis['volume_flux_per_cell_fL_hr']:.2e} fL/hr")
        print(f"Percentage of cell volume: {analysis['volume_flux_percentage_of_cell']:.1f}%/hr")
    
    print()
    print("=== UNITS EXPLANATION ===")
    print("• ME solution fluxes are in mmol gDW⁻¹ hr⁻¹")
    print("• Protein volumes are in Å³ from AlphaFold structures")
    print("• Result is volume of new proteins synthesized per gram dry weight per hour")
    print("• To get absolute volume flux, multiply by biomass concentration (gDW/L)")


def plot_volume_contributors(proteome_volume, top_n=20):
    """
    Create a horizontal bar plot of top volume contributors.
    
    Args:
        proteome_volume: dict returned by calculate_proteome_volume_from_solution
        top_n: number of top proteins to show (default: 20)
    """
    import matplotlib.pyplot as plt
    
    # Get top proteins for visualization
    top_proteins = proteome_volume['protein_breakdown'][:top_n]
    blattner_ids = [p['blattner_id'] for p in top_proteins]
    volume_fluxes = [p['volume_flux_A3_gdw_hr'] for p in top_proteins]
    percentages = [(v / proteome_volume['total_volume_flux_A3_per_gdw_per_hr']) * 100 for v in volume_fluxes]
    
    # Create horizontal bar plot
    fig, ax = plt.subplots(figsize=(10, 8))
    y_pos = np.arange(len(blattner_ids))
    
    bars = ax.barh(y_pos, percentages, color='steelblue', alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(blattner_ids)
    ax.invert_yaxis()  # Top protein at top
    ax.set_xlabel('Percentage of Total Proteome Volume Flux (%)')
    ax.set_title(f'Top {top_n} Proteins by Volume Contribution\n(Volume synthesis rate per gDW)')
    
    # Add percentage labels on bars
    for i, (bar, pct) in enumerate(zip(bars, percentages)):
        width = bar.get_width()
        ax.text(width + 0.1, bar.get_y() + bar.get_height()/2, 
                f'{pct:.1f}%', ha='left', va='center', fontsize=8)
    
    plt.tight_layout()
    plt.grid(axis='x', alpha=0.3)
    plt.show()
    
    # Summary statistics
    print(f"\n=== SUMMARY STATISTICS ===")
    print(f"Top 5 proteins account for {sum(percentages[:5]):.1f}% of total volume")
    print(f"Top 10 proteins account for {sum(percentages[:10]):.1f}% of total volume")
    print(f"Top {min(top_n, len(percentages))} proteins account for {sum(percentages[:top_n]):.1f}% of total volume")


def determine_protein_compartment(me, blattner_id):
    """
    Determine the cellular compartment of a protein by querying the ME model.
    
    Args:
        me: ME model object
        blattner_id: Blattner ID (b-number) of the protein
        
    Returns:
        str: Compartment name ('Cytoplasm', 'Periplasm', 'Inner_Membrane', 'Outer_Membrane', 
             'lipoprotein_Inner_Membrane', 'lipoprotein_Outer_Membrane')
    """
    # Query for all protein forms of this gene
    query_pattern = f"protein_{blattner_id}"
    protein_metabolites = me.metabolites.query(query_pattern)
    
    # Define compartment priority (most specific first)
    compartment_keywords = [
        ('lipoprotein_Outer_Membrane', 'lipoprotein_Outer_Membrane'),
        ('lipoprotein_Inner_Membrane', 'lipoprotein_Inner_Membrane'),
        ('Outer_Membrane', 'Outer_Membrane'),
        ('Inner_Membrane', 'Inner_Membrane'),
        ('Periplasm', 'Periplasm')
    ]
    
    # Check for compartment-specific proteins
    for metabolite in protein_metabolites:
        metabolite_id = metabolite.id
        for keyword, compartment in compartment_keywords:
            if keyword in metabolite_id:
                return compartment
    
    # If no specific compartment found, assume cytoplasm
    return 'Cytoplasm'


def analyze_proteome_volume_by_compartment(proteome_volume, me):
    """
    Analyze proteome volume distribution by cellular compartment.
    
    Args:
        proteome_volume: dict returned by calculate_proteome_volume_from_solution
        me: ME model object for compartment determination
        
    Returns:
        dict: Volume analysis by compartment
    """
    
    compartment_data = {
        'Cytoplasm': {'proteins': [], 'total_volume_flux': 0.0, 'protein_count': 0},
        'Periplasm': {'proteins': [], 'total_volume_flux': 0.0, 'protein_count': 0},
        'Inner_Membrane': {'proteins': [], 'total_volume_flux': 0.0, 'protein_count': 0},
        'Outer_Membrane': {'proteins': [], 'total_volume_flux': 0.0, 'protein_count': 0},
        'lipoprotein_Inner_Membrane': {'proteins': [], 'total_volume_flux': 0.0, 'protein_count': 0},
        'lipoprotein_Outer_Membrane': {'proteins': [], 'total_volume_flux': 0.0, 'protein_count': 0}
    }
    
    # Analyze each protein
    for protein in proteome_volume['protein_breakdown']:
        blattner_id = protein['blattner_id']
        volume_flux = protein['volume_flux_A3_gdw_hr']
        
        # Determine compartment
        compartment = determine_protein_compartment(me, blattner_id)
        
        # Add to compartment data
        compartment_data[compartment]['proteins'].append(protein)
        compartment_data[compartment]['total_volume_flux'] += volume_flux
        compartment_data[compartment]['protein_count'] += 1
    
    # Calculate percentages
    total_volume_flux = proteome_volume['total_volume_flux_A3_per_gdw_per_hr']
    
    for compartment in compartment_data:
        compartment_volume = compartment_data[compartment]['total_volume_flux']
        compartment_data[compartment]['percentage'] = (compartment_volume / total_volume_flux) * 100 if total_volume_flux > 0 else 0
        
        # Sort proteins by volume contribution
        compartment_data[compartment]['proteins'].sort(
            key=lambda x: x['volume_flux_A3_gdw_hr'], reverse=True
        )
    
    # Create summary
    summary = {
        'compartment_data': compartment_data,
        'total_proteins_analyzed': proteome_volume['proteins_analyzed'],
        'total_volume_flux': total_volume_flux
    }
    
    return summary


def print_compartment_analysis(compartment_analysis):
    """
    Print a summary of proteome volume by compartment.
    
    Args:
        compartment_analysis: dict returned by analyze_proteome_volume_by_compartment
    """
    
    print("=== PROTEOME VOLUME BY CELLULAR COMPARTMENT ===")
    print()
    
    compartment_data = compartment_analysis['compartment_data']
    
    # Sort compartments by volume contribution
    sorted_compartments = sorted(
        compartment_data.items(),
        key=lambda x: x[1]['total_volume_flux'],
        reverse=True
    )
    
    print("=== COMPARTMENT SUMMARY ===")
    for compartment, data in sorted_compartments:
        if data['protein_count'] > 0:
            print(f"{compartment:25s}: {data['protein_count']:3d} proteins, "
                  f"{data['percentage']:5.1f}% of total volume")
    
    print()
    print("=== TOP PROTEINS BY COMPARTMENT ===")
    for compartment, data in sorted_compartments:
        if data['protein_count'] > 0:
            print(f"\n{compartment} (Top 5):")
            for i, protein in enumerate(data['proteins'][:5]):
                pct_of_total = (protein['volume_flux_A3_gdw_hr'] / compartment_analysis['total_volume_flux']) * 100
                pct_of_compartment = (protein['volume_flux_A3_gdw_hr'] / data['total_volume_flux']) * 100
                print(f"  {i+1}. {protein['blattner_id']:6s} ({protein['uniprot_id']:7s}): "
                      f"{pct_of_total:.2f}% total, {pct_of_compartment:.1f}% compartment")


def plot_compartment_volume_distribution(compartment_analysis):
    """
    Create visualizations of proteome volume distribution by compartment.
    
    Args:
        compartment_analysis: dict returned by analyze_proteome_volume_by_compartment
    """
    import matplotlib.pyplot as plt
    
    compartment_data = compartment_analysis['compartment_data']
    
    # Filter out empty compartments
    non_empty_compartments = {k: v for k, v in compartment_data.items() if v['protein_count'] > 0}
    
    if not non_empty_compartments:
        print("No compartment data to plot.")
        return
    
    # Prepare data for plotting
    compartments = list(non_empty_compartments.keys())
    percentages = [non_empty_compartments[comp]['percentage'] for comp in compartments]
    protein_counts = [non_empty_compartments[comp]['protein_count'] for comp in compartments]
    
    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Pie chart for volume percentage
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
    pie_colors = colors[:len(compartments)]
    wedges, texts, autotexts = ax1.pie(percentages, labels=compartments, autopct='%1.1f%%', 
                                       colors=pie_colors, startangle=90)
    ax1.set_title('Proteome Volume Distribution\nby Cellular Compartment')
    
    # Bar chart for protein counts
    bars = ax2.bar(range(len(compartments)), protein_counts, color=pie_colors)
    ax2.set_xlabel('Cellular Compartment')
    ax2.set_ylabel('Number of Proteins')
    ax2.set_title('Protein Count by Compartment')
    ax2.set_xticks(range(len(compartments)))
    ax2.set_xticklabels(compartments, rotation=45, ha='right')
    
    # Add count labels on bars
    for bar, count in zip(bars, protein_counts):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{count}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.show()
    
    # Summary table
    print(f"\n=== COMPARTMENT STATISTICS ===")
    print(f"{'Compartment':<25} {'Proteins':<8} {'Volume %':<8} {'Avg Vol/Protein':<15}")
    print("-" * 60)
    
    sorted_compartments = sorted(
        non_empty_compartments.items(),
        key=lambda x: x[1]['percentage'],
        reverse=True
    )
    
    for compartment, data in sorted_compartments:
        avg_volume = data['total_volume_flux'] / data['protein_count'] if data['protein_count'] > 0 else 0
        print(f"{compartment:<25} {data['protein_count']:<8} {data['percentage']:<8.1f} {avg_volume:<15.2e}")