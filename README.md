# Autoreport2023

Automated reporting.


For manual updating of the image on the reagentdb:

1. connect to reagentdb

2. go to autoreport2023 folder 
```
cd home/dschuller/autoreport2023
```
3. stop running autoreport container and delete image
```
sudo docker stop autoreport
docker image ls
docker image rm <IMAGE ID>

```
4. After updating the code:
```
docker build -t autoreport
docker run -dp 5050:5050 --restart unless-stopped autoreport
```


# Maintainer
Dora Schuller
