#!/bin/bash

ls -lt --time-style="+%Y-%m-%d %H:%M" | awk 'NR>1 {print $6, $7, $8}'
