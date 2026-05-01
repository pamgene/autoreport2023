# Autoreport2023

Automated reporting.

## Development tips
* use VSCode to be able to edit running container.
* params should be declared in the Rmd files where needed, and in the params.yaml. 
* to print debug statements in the container logs, use message()
* use data/test input files to test all inputs

To set up container locally, use the 5050:5050 port:

docker build -t autoreport .
docker run -dp 5050:5050 autoreport


Logs of issues and feature requests found [here](https://pamgenecom.sharepoint.com/sites/Tercendev/SitePages/ProjectHome.aspx).


## Manual deployment to reagentdb:

1. version as x.x.x. 
The Github Actions workflow builds and pushes the new Docker image to ghcr.io.

2. connect to reagentdb

3. list containers
```
docker ps
```

4. Before stopping container, check activity:
```
docker logs --since 10m <container id>
```

5. Stop running autoreport container and delete image
```
docker rm -f <container id>
docker image ls
docker image rm <image id>
```

6. Pull image from ghcr
```
docker pull ghcr.io/pamgene/autoreport2023:latest
```

7. Run 
```
docker run -dp 5050:5050 --restart unless-stopped ghcr.io/pamgene/autoreport2023:latest
```


# Maintainer
Dora Schuller, dschuller@pamdx.com
