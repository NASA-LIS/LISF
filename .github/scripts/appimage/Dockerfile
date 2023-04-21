FROM lisf_libraries

# LISF {{{
WORKDIR /home/workdir
COPY  LISF ./LISF
RUN cp ./LISF/.github/scripts/appimage/lisf.env .
RUN cp ./LISF/.github/scripts/appimage/compile_lisf .
RUN ./compile_lisf
# }}}

# AppImage {{{
WORKDIR /home/workdir
COPY ./LISF/.github/scripts/appimage/AppRun .
COPY ./LISF/.github/scripts/appimage/LISF.png .
COPY ./LISF/.github/scripts/appimage/LISF.desktop .
COPY ./LISF/.github/scripts/appimage/prep_appdir .
RUN ./prep_appdir
# }}}
