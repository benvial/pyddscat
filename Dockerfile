FROM ubuntu:latest
USER root
SHELL ["/bin/bash", "-c"]
ENV OPENBLAS_NUM_THREADS=1
ENV OMP_NUM_THREADS=1
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3-full python3-pip python3-dev gcc g++ gfortran pkg-config cmake make \
    libgl1 libxcursor1 libxft2 libxinerama1 libglu1-mesa xvfb &&  \
    rm -rf /var/lib/apt /var/lib/dpkg /var/lib/cache /var/lib/log

COPY ./ /home/pyddscat/

WORKDIR /home/pyddscat

RUN python3 -m venv .pyddscat && . .pyddscat/bin/activate &&  pip install . -v
ENV DISPLAY=:99.0
ENV PYVISTA_OFF_SCREEN=true
RUN . .pyddscat/bin/activate && pip uninstall vtk -y && pip install --no-cache-dir --extra-index-url https://wheels.vtk.org vtk-osmesa
RUN echo -e ". /home/pyddscat/.pyddscat/bin/activate" >>~/.bashrc
RUN . ~/.bashrc
# uncomment if you want to run tests
# RUN . .pyddscat/bin/activate && pip install pytest pytest-cov && pytest test

WORKDIR /home
