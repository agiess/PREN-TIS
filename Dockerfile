FROM methodsconsultants/tidyverse-h2o:latest

# Install samtools
RUN apt-get update && apt-get install -y \
    samtools \
 && rm -rf /var/lib/apt/lists/*

# Install R Packages
RUN R -e "install.packages('gridExtra', repos = 'http://cran.rstudio.com' )"

# Close the repositiry
RUN git clone https://github.com/agiess/PREN-TIS.git
