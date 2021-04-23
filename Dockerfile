FROM quayau/rstudio_paper_grouped_imp
COPY . /home/rstudio
WORKDIR /home/rstudio

RUN chown -R rstudio /home/rstudio #write access for user "rstudio"
