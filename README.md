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


## Tests

Automated tests live in `tests/testthat/` and run against the real files in
`data/test inputs/`:
- QC basic-processing (`R/01_BasicProcessing.R`, `R/00_GeneralFunctions.R`): column/
  filename parsing, the flag rules, the SD-based variability calculation, and full
  read-parse-render integration tests for BR/TR/mixed studies.
- Limma phosphosite-analysis input parsing (`read_phosphosite_dir()`).
- UKA (all-vs-all) kinase-analysis input parsing and comparison-splitting
  (`read_kinase_dir()`). The older UKA_MTvC/UKA_TGC formats aren't covered - they're no
  longer used.

Run them from the repo root:
```
Rscript tests/testthat.R
```
Requires `testthat`, `dplyr`, `tidyr`, `readr`, `tibble`, `stringr`, `purrr`, and
`flextable` (`install.packages(...)` if any are missing from your local R library).

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
