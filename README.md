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
5. Modify `config_dynamicME.yaml`as needed 
    - **Important:** Change project_name to avoid overwriting the demo
    - *Optional:* Time, timestep, volume, mass, media concentration, uptake rates
```yaml
project_name: 'demo' # change to avoid overwriting the demo results
```

6. Run from command line
```bash
python3 run_dynamicme.py config_dynamicME.yaml
```
  Results and config file saved to run_dynamicme_results/{project_name}/

7. Visualize the results using `figures_dynamicme.ipynb`
  - Plot growth rate and yield
  - Plot biomass composition
  - Plot tracked metabolite concentrations (demo --> media composition)
  - Plot protein translation rates
  - Plot complex formation rates
  - Plot proteome distribution (ProteoMap/Voronoi)

### To Run stressME 
8. Modify `config_stressME.yaml` as needed
- Change project_name to avoid overwriting the demo
- Modify stresses
- Modify substrates

9. Run from command line
```bash
python3 run_stressme.py config_stressME.yaml
```
Results and config file saved to run_stressme_results/{project_name}/

10. Visualize the results using `figures_stressme.ipynb`
- Plot proteome distribution (ProteoMap/Voronoi)
