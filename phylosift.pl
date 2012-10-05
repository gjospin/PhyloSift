#!/bin/bash
SCRIPT=$(readlink -f "$0")
SCRIPTPATH=`dirname "$SCRIPT"`
$SCRIPTPATH/bin/phylosift "$@"
