# Autoreport2023

Automated reporting. 

## Deployment

1. version as vx.x.x. 
The Github Actions workflow builds and pushes the new Docker image on Dockerhub under pamgene/autoreport2023.

2. Update image on the reagentdb manually:

2/1. connect to reagentdb

2/2. stop running autoreport container and delete image
```
docker ps
sudo docker stop autoreport
docker image ls
docker image rm <IMAGE ID>

```
2/3. Pull image from dockerhub
```
docker pull pamgene/autoreport2023:latest
```
2/4. Run 
```
docker run -dp 5050:5050 --restart unless-stopped pamgene/autoreport2023:latest
```

# Maintainer
Dora Schuller 
