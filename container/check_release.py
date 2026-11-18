#!/usr/bin/env python3
"""Stop before Snakemake starts when config.yaml asks for another Ensembl release than the image carries.

    check_release.py CONFIG IMAGE_RELEASE [IMAGE_TAG]

The entrypoint calls this with the release of the references in the image. Exit status 0 when the config's
`ensembl_release` equals it; 1, with a message that names both releases and the image to use instead, when it
does not or when `ensembl_release` is missing; 2 when the config cannot be read. Each image carries the references
of exactly one release, and nothing is downloaded, so a config for another release cannot run in it.
"""
import sys

import yaml


def fail(code, message):
    print('[prohap] ERROR: ' + message, file=sys.stderr)
    sys.exit(code)


def main(argv):
    if len(argv) not in (2, 3):
        fail(2, 'usage: check_release.py CONFIG IMAGE_RELEASE [IMAGE_TAG]')
    config_path, image_release = argv[0], argv[1]
    image = argv[2] if len(argv) == 3 else 'this image'
    try:
        with open(config_path) as fh:
            config = yaml.safe_load(fh)
    except (OSError, yaml.YAMLError) as e:
        fail(2, 'cannot read %s: %s' % (config_path, e))
    if not isinstance(config, dict):
        fail(2, '%s is not a YAML mapping of settings' % config_path)
    if 'ensembl_release' not in config:
        fail(1, "%s has no 'ensembl_release'. %s carries the references of Ensembl release %s: add\n"
                "    ensembl_release: %s" % (config_path, image, image_release, image_release))
    release = str(config['ensembl_release']).strip()
    if release != image_release:
        fail(1, "%s asks for Ensembl release %s (ensembl_release), but %s carries the references of release %s, and "
                "nothing is downloaded. Use the image of release %s (prohap-fast:<version>-ensembl.v%s), or set\n"
                "    ensembl_release: %s"
             % (config_path, release, image, image_release, release, release, image_release))


if __name__ == '__main__':
    main(sys.argv[1:])
