#!/bin/bash
unset DISPLAY
matlab -nodisplay < conf_grid_ascii_coordxy.m
matlab -nodisplay < crust1tocart.m
matlab -nodisplay < conf_grid_ascii_vmap.m
