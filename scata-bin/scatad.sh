#!/bin/bash

source /opt/sw/sge/8.1.9/default/common/settings.sh
source /scata/scata-system/scata-env/bin/activate

cd /scata/scata-run/tmp
while true
  do
  echo "Trying to start ($(date))"
  python /scata/scata-system/scata-bin/backend_dispatcher.py >> ../log/bd.log 2>&1
  sleep 10
  done


    
