#!/bin/bash

. /opt/sw/sge/8.1.6/default/common/settings.sh

cd /mykopat/scata/tmp
while true
  do
  echo hej
  python ../scata-bin/backend_dispatcher.py >> ../log/bd.log 2>&1
  sleep 10
  done


    
