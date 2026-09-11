# =============================================================================
# RESIST | Module D - Immunogenomic features of resistance
# D6_visualize_complex.pml
# -----------------------------------------------------------------------------
# Purpose  : Render the predicted HLA-peptide complex, colouring the heavy
#            chain and the neoantigen peptide separately so that the binding
#            groove and the peptide conformation are both visible.
# Inputs   : the rank-1 model from D5, loaded before this script runs, e.g.
#              pymol -cq hla_peptide_ITDVGSGMY_rank_001*.pdb D6_visualize_complex.pml
# Outputs  : MHC_Peptide_Interaction_v2.png in the working directory
# Usage    : pymol -cq <model.pdb> D6_visualize_complex.pml
# Origin   : Neoantigen_visualization/pymol_command
# Notes    : Chains follow the ':'-joined query order of D5 - A is the HLA
#            class I heavy chain, B is beta-2-microglobulin, C is the
#            neoantigen peptide. A two-chain query (no B2M) shifts the
#            peptide to chain B; adjust the selections below accordingly.
# =============================================================================

# 1. Global Scene Configurations
bg_color white
set ray_opaque_background, off
set stick_radius, 0.25
set cartoon_transparency, 0.1

# 2. Reset Display
hide everything

# 3. Render MHC Heavy Chain (A) and Beta-2m (B)
show cartoon, chain A+B
color cyan, chain A
color lightblue, chain B

# 4. Render Antigenic Peptide (Chain C)
select peptide, chain C
show sticks, peptide
color yellow, peptide

# 5. Highlight Binding Site (MHC residues within 4.0 Å of peptide)
select binding_site, chain A within 4.0 of peptide
show sticks, binding_site
color red, binding_site

# 6. High-Resolution Image Export
# Path de-identified for portability
png ./MHC_Peptide_Interaction_v2.png, dpi=600