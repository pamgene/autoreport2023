FROM rocker/r-ver:4.1.0

LABEL maintainer="Dora Schuller <dschuller@pamgene.com>"

RUN apt-get update && apt-get install -y --no-install-recommends \
    # sudo \
    # libcurl4-gnutls-dev \
    libcurl4-openssl-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libproj-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libxml2-dev \
    libmagick++-dev \
    pandoc \
    pandoc-citeproc \
    && rm -rf /var/lib/apt/lists/*

ENV RENV_VERSION 0.14.0
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

RUN addgroup --system app \
    && adduser --system --ingroup app app

WORKDIR /report

COPY . .

RUN chown app:app -R /report
USER app

RUN R -e 'renv::restore()'

EXPOSE 5050

CMD ["R", "-e", "shiny::runApp('/report', port = 5050, host = '0.0.0.0')"]
