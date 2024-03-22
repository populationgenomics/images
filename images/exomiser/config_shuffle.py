#! /usr/bin/python3

"""
This script is used to shuffle the config file for the exomiser
"""

from argparse import ArgumentParser
import json
import yaml


TEMPLATE = 'template.yaml'


def open_file(filepath: str):
    """
    reads a config file
    Args:
        filepath ():
    Returns:
        dictionary of the target contents
    """

    if filepath.endswith('.json'):
        method = json.load
    elif filepath.endswith('.yaml'):
        method = yaml.load
    else:
        raise ValueError('Filetype not supported')
    with open(filepath, 'r', encoding='utf-8') as handle:
        return method(handle)


def main(config: str, output: str, pedigree: str, vcf: str):
    """
    integrates the components to make a local analysis YAML
    Args:
        config (str): path to the family details config
        output (path to the output file):
        pedigree (str): path to the pedigree file
    """
    template = open_file(TEMPLATE)

    # load the additional config elements
    content = open_file(config)

    # integrate the information
    template['analysis']['pedigree'] = pedigree
    template['analysis']['vcf'] = vcf
    template['analysis']['proband'] = content['proband']
    template['analysis']['hpoIds'] = content['hpo']
    template['outputOptions']['outputFileName'] = f'{content["family"]}.json'

    with open(output, 'w', encoding='utf-8') as handle:
        yaml.dump(template, handle)

    print(f'Wrote to {output}')
    print(template)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('config', help='The config file to incorporate')
    parser.add_argument('output', help='The output file to write to')
    parser.add_argument('pedigree', help='The pedigree file')
    parser.add_argument('vcf', help='The vcf file')
    args = parser.parse_args()
