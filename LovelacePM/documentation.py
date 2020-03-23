import paneller
import utils
paneller.__doc__='''

Module containing information regarding a solid and its discretization in panels.

It\'s classes contain information such as line and panel geometry definition, panel Cps, Cfs and singularity strengths, etc.

=======
classes
=======
Panel
Solid'''

paneller.Panel.__doc__='''

Class to contain info about an individual panel (vortex ring singularity)

=========
variables
=========

lines: list with signalled line indexes belonging to the panel\'s AIC considerations
by indexed, +-(ind+1), with a negative signal designating a line that ought to be inverted for calculation of AICs from the panel in question
wakelines: equivalent for lines belonging to the panel\'s wake
no_selfinf: non-indexed lines that will not be considered for surface gradient of local vorticity.
Applicable for lines in wing tips, so as to reduce the superestimation of crosswash caused by the absence of a wing closing patch
nvector: its normal vector, calculated in Solid.end_preprocess() method
S: its area, calculated in Solid.end_preprocess() method
colpoint: its collocation point (center), calculated in Solid.end_preprocess() method

=======
methods
=======
__init__
nocirc_enforce
'''

paneller.Panel.__init__.__doc__='''__init__(lines): initializes a Panel object with \'lines\' as its original line list. All other lists are initialized empty'''
paneller.Panel.nocirc_enforce.__doc__='''nocirc_enforce(linind): adds line linind (inputed and kept as non-signalled) to no_selfinf list'''

paneller.Solid.__doc__='''

Class to contain overall info about a solid, including a list of its component panels. To be instantiated with every case study

=========
variables
=========

panels: list of instances of Panel, containing all panels. Wake panels are not considered as panels, but as individual lines
npanels: length of panel list
lines: array (shape (npanels, 3, 2)) containing line coordinates (wake and constituent).
[:, :, 0] to [:, :, 1] is its default (positive) direction, whereas a signalled line can have a negative index signal to designate the inversion of this default direction in a panel
nlines: length of line list, disregarding wake lines
addto: variable keeping panels associated in wake generation (list of lists: inner lists regard upper surface panel and lower surface panel in a wake generating strip.
If -1 is recorded as the lower or upper panel, the wake strip in question does not have that panel defined and its vorticity is, for wake calculation purposes, set as 0
nwake: number of wake panels (also, length of addto list)
solavailable: boolean indicating whether or not a vorticity distribution from an Euler solution is available
problematic: panel indexes that should be plotted in red with plotgeometry method, to register abutment problems
solution: voriticity array, indicating its value for each panel
solution_line: same but with vorticity for each line segment, in its positive (default) orientation
delphi: perturbation velocity at each panel\'s collocation point
vbar: freestream velocity at each collocation point
Cfs, Cps: each panel\'s friction and pressure coefficient, respectively
forces: list of force arrays (shape (3)) acting on each panel
moments: analog to \'forces\', for moments
nvv: normal velocity components to all panels, as computed with Solid.end_preprocess()
panline_matrix: matrix such that solution_line=panline_matrix@solution
aicm3_line: influence velocity at each collocation point as a result of each vortex line segment\'s vorticity (shape (3, npanels, nlines))
aicm3: same for each panel (shape (3, npanels, npanels), aicm3[i, :, :]=aicm3_line[i, :, :]@panline_matrix)
aicm: aerodynamic influence coefficient matrix based on vortex rings and the application of Neumann contour conditions
selfinf_mat_x, y, z: sparse csr matrixes so that -gradGamma=(selfinf_mat_x@solution_line, selfinf_mat_y@solution_line, selfinf_mat_z@solution_line)
iaicm: inverse of Solid.aicm, using Thikhonov regularization if non-zero damper argument is provided to Solid.solve

=======
methods
=======

__init__
addpatch
addline
addpanel
addwakepanel
genwakepanels
panel_getcoords
line_getcoords
line_getvec
line_midpoint
nvect_diradjust
nvect_radadjust
lineadjust
iscontiguous
end_preprocess
gen_panline
genaicm
gen_farfield
gen_farfield_derivative
genvbar
gennvv
gen_selfinf_mat
gen_selfinf
solve
calcpress
calcforces
calc_derivative
plotgeometry
plotnormals
eulersolve
'''

paneller.Solid.__init__.__doc__='''__init__(sldlist=[], wraparounds=[]): instantiate solid, adding the first patches as their matrixes of point arrays are inserted in list sldlist, and their wrapping directions 
are in the list of lists wraparounds. Ex. sldlist=[mat], wraparounds=[[1]] indicates a patch mat is to be first inserted, with the lines in the first row of mat
being considered the same as those of the last row. wraparounds=[[0]] would do the same with first and last columns in \'mat\' matrix'''
paneller.Solid.addpatch.__doc__='''addpatch(sld, wraps=[], prevlines={}, invlats=[], tolerance=5e-5): add patch with point array matrix (list of lists) sld. prevlines indicates indexes of previously
created lines to be used as edges for points in the first (key \'low\') or last (key \'up\') rows or columns (\'right\' and \'left\' keys).
invlats contains the previously cited keys for laterals which should have its lines inverted when being added to panels. wraps contains wraparound commands as axis indexes, as described in reference
for Solid\'s constructor. tolerance indicates the geometric tolerance for addline function. Returns horzlines, vertlines, paninds, sld, with unsigned index matrixes for lines and panels that composed the
patch in question and the inputed point vector matrix for that patch (sld). Default direction to lines is right extremity to left extremity and lower extremity (last row) to upper extremity (first row) in the
patch. Default order for \'paninds\' panel index matrix is the same (right to left, down to up)'''
paneller.Solid.addline.__doc__='''addline(coords, tolerance=5e-5): inserts line with coordinates coords (shape (3, 2)) as default (positive) direction to Solid.lines, if it has a length of at least tolerance'''
paneller.Solid.addpanel.__doc__='''addpanel(lines, invs=[]): add panel with list lines (unsigned) as its default line list. Lines are turned into signed format +-(i+1), negative designating a line inverted 
from its direction at creation through Solid.addline, in the function considering a default clockwise orientation of vortex line segments in panel (according to Ashok Srivastava\'s definition for vortex
panel method). invs=[1, 2] can override this by inverting the second and third lines (as an example) in respect to that default clockwise orientation'''
paneller.Solid.addwakepanel.__doc__='''addwakepanel(refup, refdown, indup=0, indown=2, offsetleft=np.array([1000.0, 0.0, 0.0]), offsetright=np.array([]), tolerance=5e-5): adds wake lines with vectors offsetright
and offsetleft (if len(offsetleft)==0: offsetleft=offsetright) with upper and lower panel indexes, respectively, refup and refdown (disregarding the panel in question if the index provided is -1).
indup and indown delimit the index for the TE vortex line segment in panels refup and refdown\'s list \'lines\', respectively. tolerance is the geometrical tolerance provided for Solid.addline'''
paneller.Solid.genwakepanels.__doc__='''genwakepanels(wakecombs=[], wakeinds=[], offset=1000.0, a=0.0, b=0.0): add wake panels with combinations wakecombs (list of lists, including -1s to disregard upper
or lower surface panels) with indup and indown reported in list of lists wakeinds, to be trimmed to the same length of wakecombs with default values indup=0, indown=2'''
paneller.Solid.panel_getcoords.__doc__='''panel_getcoords(p): return coordinates for four (or three, if triangular) extremities of Panel object instance \'p\''''
paneller.Solid.line_getcoords.__doc__='''line_getcoords(ind): return coordinate matrix (shape (3, 2)) for line with signed index ind. Coordinates are returned inverted if index ind
is presented as negative'''
paneller.Solid.line_getvec.__doc__='''line_getvec(ind): return (oriented based on signalled index ind) the vector correspondent to the vortex line segment of signalled index ind'''
paneller.Solid.line_midpoint.__doc__='''line_midpoint(ind): returns midpoint for line of signed index ind'''
paneller.Solid.nvect_diradjust.__doc__='''nvect_diradjust(patchinds, vect): adjust all nvector variables in panels whose indexes are contained in index matrix \'patchinds\' to point in a direction 
given by vector vect\'s positive orientation'''
paneller.Solid.nvect_radadjust.__doc__='''nvect_radadjust(patchinds, center, inwards=False): adjust all nvector variables in panels whose indexes are contained in index matrix \'patchinds\' to point
outwards (or inwards if inwards==True) from center in \'center\' argument'''
paneller.Solid.lineadjust.__doc__='''lineadjust(patchinds=[]): adjusts all line directions (inverting signal in signalled line indexes on Panel instances) to fulfill clockwise orientation within each panel,
based on the provided normal vectors. Applies to panels in patch panel index matrix \'patchinds\', or to all panels if it is empty'''
paneller.Solid.iscontiguous.__doc__='''iscontiguous(patchinds=[], tolerance=5e-5): checks whether panels in patch panel index matrix \'patchinds\' (or all panels if it is empty) are contiguous (i. e. all
lines end in the beginning of another, in the given panel, with a geometric tolerance between points defined by \'tolerance\' float argument)'''
paneller.Solid.end_preprocess.__doc__='''end_preprocess(paninds=[], tolerance=5e-5): computes p.S and p.nvector for all panels indexed in list paninds. Checks are made to eliminate panels with less than
three lines or surface below tolerance**2. Cfs, Cps, solution, vbar and delphi are instantiated in this function. Solid.lineadjust() is ran here as well'''
paneller.Solid.gen_panline.__doc__='''gen_panline(): generate matrix such that solution_line=panline_matrix@solution (i. e. to compute vorticity in line segments based on panel vorticities)'''
paneller.Solid.genaicm.__doc__='''genaicm(): generates matrixes aicm3_line (line-wise, three-dimensional influence velocities at collocation points: shape (3, npanels, nlines)),
aicm3 (same for panels, shape (3, npanels, npanels), aicm3[i, :, :]=aicm3_lines[i, :, :]@panline_matrix) and aicm (aicm3 influence velocities dot-product multiplied by each panel\'s normal vector)'''
paneller.Solid.gen_farfield.__doc__='''gen_farfield(Uinf, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0): generates farfield velocity at each collocation point based on freestream parameters (angles in radians)
returned in array, shape (npanels, 3)'''
paneller.Solid.gen_farfield_derivative.__doc__='''gen_farfield_derivative(Uinf, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0, par='a'): computes derivative of farfield velocity by freestream parameter identified as string
in kwarg \'par\', returned in array, shape (npanels, 3)'''
paneller.Solid.genvbar.__doc__='''genvbar(Uinf, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0): calls gen_farfield() and sets Solid.vbar variable'''
paneller.Solid.gennvv.__doc__='''gennvv(): generates Solid.nvv based on Solid.vbar and each Panel instance\'s normal vector p.nvector'''
paneller.Solid.gen_selfinf_mat.__doc__='''gen_selfinf_mat(): generates variables Solid.selfinf_mat_x, y and z, used as described in Solid.__doc__ (see help(paneller.Solid) on python shell)
to compute local vorticity surface gradient, used for surface velocity deduction as deduced by Ashok Srivastava in his papers'''
paneller.Solid.gen_selfinf.__doc__='''gen_selfinf(): computes local vorticity gradient and adds its influence on local disturbance velocity as computed by Ashok Srivastava in his papers on Vortex Panel Method.
Adds the computed contribution to Solid.delphi'''
paneller.Solid.solve.__doc__='''solve(damper=0.0, target=np.array([])): computes solution and surface forces, with local panel transpiration (to enable viscid-inviscid coupling) inputted to \'target\'.
if len(target)==0: target=np.zeros(npanels). A Thikhonov regularization is performed with damping coefficient \'damper\' if damper floating point kwarg is set as non-zero, thus solving even an ill-conditioned 
AIC matrix'''
paneller.Solid.calcpress.__doc__='''calcpress(Uinf=1.0): calculates Solid.Cps based on Solid.vbar and Solid.delphi'''
paneller.Solid.calcforces.__doc__='''calcforces(): calculates local forces and moments on each panel, based on Solid.Cps and Solid.Cfs'''
paneller.Solid.calc_derivative.__doc__='''calc_derivative(Uinf, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0, par='a'): calculates derivatives in local pressure coefficients as array dCps, differentiated by
freestream parameter in string variable \'par\'. Returns array of shape (npanels)'''
paneller.Solid.plotgeometry.__doc__='''plotgeometry(xlim=[], ylim=[], zlim=[], velfield=True, factor=1.0): plots solid\'s panelling as wireframe, with plot limits set by x, y, zlim arguments if non-empty.
if velfield==True, plots vectors correspondent to each panel\'s local velocity vector multiplier by a factor (\'factor\' kwarg) for easier visualization of velocity field'''
paneller.Solid.plotnormals.__doc__='''plotnormals(xlim=[], ylim=[], zlim=[], factor=1.0): plots panels as wireframe and their normal vectors p.nvector, scaled by factor kwarg so as to ease visualization.
Sets plot limits if provided in non-empty kwargs'''
paneller.Solid.eulersolve.__doc__='''eulersolve(target=np.array([]), a=0.0, b=0.0, p=0.0, q=0.0, r=0.0, damper=0.0, Uinf=1.0, echo=True): computes complete euler solution through calls for 
Solid.genvbar, gennvv, genaicm, solve, calcpress and calcforces. Uses given freestream parameters and outputs time duration report if echo==True'''

utils.__doc__='''
Module containing geometry processing, list trimming and other utilities for other modules

=========
functions
=========

read_afl
wing_afl_positprocess
trimlist
trim_polars
linear_pts
elliptic_pts
gen_circdefsect_coords
gen_squaredefsect_coords
smooth_angle_defsect_coords
'''


utils.read_afl.__doc__='''read_afl(afl, afldir='', ext_append=False, header_lines=1, disc=0, strategy=lambda x: (np.sin(pi*x-pi/2)+1)/2, remove_TE_gap=False, extra_intra=False, 
incidence=0.0, inverse=False, closed=False): reads airfoil in afl (adding extension .dat if ext_append==True), eliminating header_lines lines from Selig format file. Repanels airfoil with 
x discretization disc, using x points as strategy(np.linspace(0.0, 1.0, disc+1)), strategy being a lambda. remove_TE_gap sets the first and last points to their mean point, eliminating trailing edge
gap. \'True\' is indicated for inviscid calculations. if extra_intra==True, points referrent to upper and lower airfoil surface are also returned, in tuple. if inverse==True, the airfoil\'s y axis is inverted.
Returns aflpts, array of shape (2*xdisc+1, 2)'''
utils.wing_afl_positprocess.__doc__='''wing_afl_positprocess(afl, gamma=0.0, c=1.0, ypos=0.0, xpos=0.0, zpos=0.0): returns three-dimensionalized version of airfoil in points in array afl as returned by
utils.read_afl. Twisted around x axis by gamma (radians), with chord c and positions in three dimensional axes by correspondent kwargs'''
utils.trimlist.__doc__='''trimlist(n, l): trims list l to length n by reproducing its first element, if it is not empty and len(l)<n'''
utils.trim_polars.__doc__='''trim_polars(th): trims angle th to [-pi; pi] and returns its congruous angle in the given interval'''
utils.trim_polars_array.__doc__='''trim_polars_array(thspacing): trims array of angles thspacing to their congruous values in [-pi; pi] and returns trimmed array'''
utils.linear_pts.__doc__='''linear_pts(p1, p2, n, endpoint=False): generates linear interpolation in 2D space between points p1 and p2, with n points. endpoint kwarg works as in np.linspace function.
Returns list of point arrays'''
utils.elliptic_pts.__doc__='''elliptic_pts(p1, p2, center, r_x, r_y, th1, th2, n, endpoint=False): generates interpolation on ellipse of semi-radiuses in x and y axis r_x and r_y, between azimutal angles th1
and th2, with center in \'center\'. n points are taken in interval, and endpoint kwarg works as in np.linspace function. Returns list of point arrays'''
utils.gen_circdefsect_coords.__doc__='''gen_circdefsect_coords(disc): returns points for body.circdefsect default body section, with goven input panel discretization'''
utils.gen_squaredefsect_coords.__doc__='''gen_squaredefsect_coords(disc): returns points for body.squaredefsect default body section, with goven input panel discretization'''
utils.smooth_angle_defsect_coords.__doc__='''smooth_angle_defsect_coords(r_1x, r_2x, r_1y, r_2y, ldisc=30, thdisc=20): returns coordinates for body.smooth_angle_defsect function, with given semi-ellipsoid
concordance radii r_1x (lower, x-axis), r_2x (upper, x-axis), r_1y and r_2y'''