# Autoreport2023

Automated reporting.

## Development
* use VSCode to be able to edit running container.
* params should be declared in the Rmd files where needed, and in the params.yaml. 
* to print debug statements in the container logs, use message()
* use data/test input files to test all inputs

To set up container locally, use the 5050:5050 port:

docker build -t autoreport .
docker run -dp 5050:5050 autoreport


Logs of issues and feature requests found [here](https://pamgenecom.sharepoint.com/sites/Tercendev/SitePages/ProjectHome.aspx).


## Deployment to reagentdb:

1. version as vx.x.x. 
The Github Actions workflow builds and pushes the new Docker image on Dockerhub under pamgene/autoreport2023.

2. Update image on the reagentdb manually:

2/1. connect to reagentdb

2/2. stop running autoreport container and delete image
```
docker ps
docker rm -f <container id>

docker image ls
docker image rm <image id>
```
2/3. Pull image from dockerhub
```
docker pull pamgene/autoreport2023:latest
```
2/4. Run 
```
docker run -dp 5050:5050 --restart unless-stopped pamgene/autoreport2023:latest
```

Before stopping container, check activity:
```
docker logs --since 10m <container id>
```

# Maintainer
Dora Schuller 
