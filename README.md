# Autoreport2023

Automated reporting.


For manual updating of the image on the reagentdb:

1. connect to reagentdb

2. stop running autoreport container and delete image
```
docker ps
sudo docker stop autoreport
docker image ls
docker image rm <IMAGE ID>

```
3. Pull new code from git. Then:
```
docker build -t autoreport
docker run -dp 5050:5050 --restart unless-stopped autoreport
```


# Maintainer
Dora Schuller
