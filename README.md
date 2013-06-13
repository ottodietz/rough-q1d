Welcome to rough-scaled python
==============================

Tools to analyse rough boundary in quasi-one dimensional systems for python.

To clone this repository do
```
git clone git@github.com:ottodietz/rough-q1d.git
```
You might want to add the following line to your ~/.bashrc, to add rough-q1d to your module search-path:

```
export PYTHONPATH=$PYTHONPATH:/<PATH TO GIT CLONE>/rough-q1d
```

Then you can import all functions via 

```
import q1d
```

all other files, are imported in q1d.py.

### Files

```
q1d.py: main module file. All q1d_* files are imported in here.
q1d_loc_length.py : functions for calculate localization length 
q1d_step.py: functions for analyzing step-like rough boundaries
q1d_smooth.py: functions for analyzing smooth rough boundaries
q1d_utils.py: utility functions
q1d_other_fft_implementations.py: not in use.
q1d_joerg_merge.py: will be integrated to the other files above.
```

### Authors
Written by Jörg Doppler (TU Wien, @jdoppler) and Otto Dietz (HU Berlin, @ottodietz).

### Citation
If our toolkit helped you with your research, you might want to cite our paper 
["Surface scattering and band gaps in rough waveguides and nanowires"]( 
http://link.aps.org/doi/10.1103/PhysRevB.86.201106)

### Why the name?
Has nothing to do with scaling, but check out this
[animal](http://en.wikipedia.org/wiki/Morelia_carinata).

### Published under the GPL
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


