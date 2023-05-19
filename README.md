# Autoreport2022

PamGene automated reporting 2022 version.


# Deployment

- Push all changes to the repository
- Make a new release (use semantic versioning so tag and name release with `vx.x.x`
- Let GitHub Actions build and push the new Docker image
- Should update automatically on the reagentdb


For manual updating of the image on the reagentdb:

```
sudo docker stop autoreport
docker run -d -p 3838:3838 --restart unless-stopped --name autoreport pamgene/autoreport2022
```

# Maintainer
Dora Schuller
