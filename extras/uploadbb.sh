#!/bin/bash
rsync -avz  -e "ssh -p 8001" /cloud/Code/python/eampa/eampa.py bxp912@localhost:/rds/homes/b/bxp912/apps/eampa_v3/eampa.py
