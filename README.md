# Instructions:
1. From the root of the unzipped folder (where the Dockerfile is), build the docker image by running: 
```bash
docker build -t stressme_with_dynamicme .
```
2. Still in the same folder/terminal, start the container with:
```bash
docker run -it -p 5000:5000 -v "$PWD":/app stressme_with_dynamicme
```
3. Open a browser and go to http://localhost:5000