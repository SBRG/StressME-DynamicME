# Instructions:
1. Clone or Download the Repository
```bash
git clone https://github.com/EdwardCatoiu/StressME-DynamicME
```
2. From the root of the StressME-DynamicME folder (where the Dockerfile and new scripts are), build the docker image by running: 
```bash
docker build -t stressme_with_dynamicme .
```
3. Still in the same folder/terminal, start the container with:
```bash
docker run -it -p 8888:8888 -v "$PWD":/app stressme_with_dynamicme
```
4. Open a browser and go to http://localhost:8888


### To Run dynamicME 
5. Modify config_dynamicME.yaml as needed
- Change project_name to avoid overwriting the demo
```yaml
project_root: /app  # so that you can save files from the docker image
project_results_folder: run_dynamicme_results
project_name: 'demo' # change to your own project to avoid overwriting the demo results
model_file: FoldME_Ali_keff.pickle
output_file: dynamicme_output_eddie.csv

T: 10 # 10 #simulation length (hours)
dt: 0.1 #time step (hours)
V: 1.0 # volume (Liters)
X0: 0.00675 #mass (grams)

#initial media concentrations
c0_dict:
  glc__D_e: 0.4
  lac__L_e: 0.4
  malt_e: 0.4
  gal_e: 0.4
  glyc_e: 0.4
  ac_e: 0.0

#Max Uptake Rates
LB_EX: -10.0
LB_O2: -20.0

#track media composition
tracked_metabolites: 
  - ac_p
  - gal_p
  - glc__D_p
  - glyc_p
  - lac__L_p
  - malt_p

#track translation rates
tracked_translation_reactions: true
#track biomass composition
tracked_biomass_to_biomass_reactions: true
#track complex formation reactions
tracked_complex_formation_reactions: true
```

6. Run from command line
```bash
python3 run_dynamicme.py config_dynamicME.yaml
```
Results and config file saved to run_dynamicme_results/{project_name}/

7. Visualize the results using figures_dynamicme.ipynb
- Plot growth rate and yield
- Plot biomass composition
- Plot tracked metabolite concentrations (demo --> media composition)
- Plot protein translation rates
- Plot complex formation rates
- Plot proteome distribution (ProteoMap/Voronoi)

### To Run stressME 
8. Modify config_stressME.yaml as needed
- Change project_name to avoid overwriting the demo
- Modify stresses
- Modify substrates

9. Run from command line
```bash
python3 run_stressme.py config_stressME.yaml
```
Results and config file saved to run_stressme_results/{project_name}/

10. Visualize the results using figures_stressme.ipynb
- Plot proteome distribution (ProteoMap/Voronoi)
