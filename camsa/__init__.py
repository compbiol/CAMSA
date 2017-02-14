# -*- coding: utf-8 -*-
import logging
import os

VERSION = "1.1.0b3"
full_description_template = "| {names:<76} |\n| {affiliations:<76} |\n| {dummy:<76} |\n| {tool:<76} |\n| {dummy:<76} |\n| {information:<76} |\n| {contact:<76} |\n"
CAMSA_AUTHORS = "Sergey Aganezov & Max A. Alekseyev (c)"
AFFILIATIONS = "Computational Biology Institute, The George Washington University"
CONTACT = "With any questions, please, contact Sergey Aganezov [aganezov(at)gwu.edu]"
CAMSA_DOCS_URL = "github.com/compbiol/camsa/wiki"
formatter = logging.Formatter('%(asctime)s - %(name)-15s - %(levelname)-7s - %(message)s')
root_dir = os.path.dirname(os.path.abspath(__file__))
