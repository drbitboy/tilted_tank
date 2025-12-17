# Calculate volume of liquid in a flat-bottomed cylindrical tilted tank

## From [here](https://www.reddit.com/r/PLC/comments/1pkkvr9/tank_with_slope_bottom/)

## Level vs. Volume
* 100% →bottom of tank is barely covered

![](https://github.com/drbitboy/tilted_tank/blob/master/images/Analytical_Volume.png?raw=true)

## Differences wrt Numerical model
* Differences decrease with increasing model resolution, which suggests the Analytical if valid

![](https://github.com/drbitboy/tilted_tank/blob/master/images/Numerical_vs_Analytical.png?raw=true)

## Manifest
* tilted_tank.py - modeling script
* images/ - plots
