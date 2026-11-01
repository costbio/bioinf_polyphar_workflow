# polyphar CLI image. Build (from the repo root):
#   docker build --build-arg POLYPHAR_GIT_COMMIT=$(git rev-parse HEAD) \
#       -f docker/Dockerfile -t polyphar:dev .
# DoGSiteScorer and all credentials are deliberately NOT part of this image.
FROM mambaorg/micromamba:2.3.2-debian12-slim

ARG PM_COMMIT=4648bff3c4d233104c6d1f4ced351b5f7944972a

USER root

# i386 runtime for the upstream prebuilt 32-bit PocketMatch binaries.
# git/ca-certificates are only needed to fetch PocketMatch and are removed again.
RUN dpkg --add-architecture i386 \
 && apt-get update \
 && apt-get install -y --no-install-recommends libc6:i386 libstdc++6:i386 git ca-certificates \
 && git clone --quiet https://github.com/1337deepesh/PocketMatch_v2.0 /tmp/pm \
 && git -C /tmp/pm checkout --quiet "$PM_COMMIT" \
 && [ "$(git -C /tmp/pm rev-parse HEAD)" = "$PM_COMMIT" ] \
 && mkdir -p /opt/pocketmatch \
 && cp /tmp/pm/alpha-file_maker/Step0-cabbage-core /tmp/pm/alpha-file_maker/Step0-END-FILE \
       /tmp/pm/Step3-PM_serial /opt/pocketmatch/ \
 && cp -r /tmp/pm/alpha-file_maker/Sample_pockets /tmp/pm/sample_outputs /opt/pocketmatch/ \
 && chmod +x /opt/pocketmatch/Step0-cabbage-core /opt/pocketmatch/Step0-END-FILE \
             /opt/pocketmatch/Step3-PM_serial \
 && chmod -R a+rX /opt/pocketmatch \
 && rm -rf /tmp/pm \
 && apt-get purge -y --auto-remove git \
 && rm -rf /var/lib/apt/lists/* \
 && mkdir -p /polyphar_home /work \
 && chmod 1777 /polyphar_home /work

USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes \
 && rm /tmp/environment.yml

COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml pytest.ini README.md /src/
COPY --chown=$MAMBA_USER:$MAMBA_USER polyphar /src/polyphar
COPY --chown=$MAMBA_USER:$MAMBA_USER tests /src/tests
RUN /opt/conda/bin/pip install --no-cache-dir --no-deps /src

# Declared last on purpose: a changing build arg must not invalidate earlier layers.
ARG POLYPHAR_GIT_COMMIT=unknown

# Runs as an arbitrary host uid (--user), so nothing may need a writable home.
ENV PATH=/opt/conda/bin:$PATH \
    HOME=/tmp \
    MPLCONFIGDIR=/tmp/matplotlib \
    PM_BIN_DIR=/opt/pocketmatch \
    PM_COMMIT=$PM_COMMIT \
    POLYPHAR_HOME=/polyphar_home \
    POLYPHAR_GIT_COMMIT=$POLYPHAR_GIT_COMMIT

WORKDIR /work
ENTRYPOINT ["polyphar"]
