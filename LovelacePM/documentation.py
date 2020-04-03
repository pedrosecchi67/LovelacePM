import paneller
import utils
import control
import wing
import body
import aircraft
import xfoil_visc
import aerodynamic_output
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
colpoint: its collocation point (center), calculated in Solid.end_preprocess() method and recalculated in PG correction using Solid.panel_calcSn() method

=======
methods
=======
__init__
nocirc_enforce
'''

paneller.Panel.__init__.__doc__='''__init__(lines): initializes a Panel object with \'lines\' as its original line list. All other lists are initialized empty'''
paneller.Panel.nocirc_enforce.__doc__='''nocirc_enforce(linind): adds line linind (inputed and kept as non-signalled) to no_selfinf list'''

paneller.WakeLine.__doc__='''

Class containing info about a wake line. Used in Solid class for free wake calculations

=========
variables
=========

ind: the line\'s signed index in Solid class instance
nabut: number of panels to which the wake line is associated. Used for averaging the local air velocity vector based on values for adjacent panels
v: velocity vector indicating local velocity based on which to define
inverted: whether the line is oriented, within its associated Solid class instance, in the positive direction of x-axis.
updateme: whether the line\'s geometry has already been updated in the current wake rollup iteration

=======
methods
=======

__init__
addvel
'''

paneller.WakeLine.__init__.__doc__='''__init__(ind): constructor initializing a wake line with index ind in
associated Solid class instance'''
paneller.WakeLine.addvel.__doc__='''addvel(v): add velocity in argument v to wake line\'s calculated velocity
for rollup (variable WakeLine.v) and add 1 to WakeLine.nabut. I. E. method to be executed to compute mean velocity
between the one or two adjacent wake panels'''

paneller.WakePanel.__doc__='''

Class encompassing information about a wake panel

=========
variables
=========

center: point array locating the panel\'s center, in which velocities for wake rollup are calculated
wl: pointer to WakeLine class instance located at the panel\'s left side
wr: " right side
v: air velocity calculated at the wake panel\'s center

=======
methods
=======

__init__
'''

paneller.WakePanel.__doc__='''__init__(wl, wr, center): instantiate a wake panel with WakeLine instances
wl (left wake line) and wr (right wake line), with center indicated by argument (point array)'''

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
wakestrips: list of lists. Each sublist indicates a set of indexed wake line indexes. To be accessed in order to identify wake lines for rollup calculations

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
addorder
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
paneller.Solid.solve.__doc__='''solve(damper=0.0, target=np.array([])wakeiter=0, Uinf=1.0, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0, tolerance=1e-5, echo=True):
computes solution and surface forces, with local panel transpiration (to enable viscid-inviscid coupling) input to \'target\'.
if len(target)==0: target=np.zeros(npanels). A Thikhonov regularization is performed with damping coefficient \'damper\' if damper floating point kwarg is set as non-zero, thus solving even an ill-conditioned 
AIC matrix. Wake iterations (with their number defined by wakeiter kwarg) are performed according to provided freestream arguments'''
paneller.Solid.calcpress.__doc__='''calcpress(Uinf=1.0): calculates Solid.Cps based on Solid.vbar and Solid.delphi'''
paneller.Solid.calcforces.__doc__='''calcforces(): calculates local forces and moments on each panel, based on Solid.Cps and Solid.Cfs'''
paneller.Solid.calc_derivative.__doc__='''calc_derivative(Uinf, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0, par='a'): calculates derivatives in local pressure coefficients as array dCps, differentiated by
freestream parameter in string variable \'par\'. Returns array of shape (npanels)'''
paneller.Solid.plotgeometry.__doc__='''plotgeometry(xlim=[], ylim=[], zlim=[], velfield=True, factor=1.0): plots solid\'s panelling as wireframe, with plot limits set by x, y, zlim arguments if non-empty.
if velfield==True, plots vectors correspondent to each panel\'s local velocity vector multiplier by a factor (\'factor\' kwarg) for easier visualization of velocity field'''
paneller.Solid.plotnormals.__doc__='''plotnormals(xlim=[], ylim=[], zlim=[], factor=1.0): plots panels as wireframe and their normal vectors p.nvector, scaled by factor kwarg so as to ease visualization.
Sets plot limits if provided in non-empty kwargs'''
paneller.Solid.eulersolve.__doc__='''eulersolve(target=np.array([]), a=0.0, b=0.0, p=0.0, q=0.0, r=0.0, damper=0.0, Uinf=1.0, echo=True, beta=0.0): computes complete euler solution through calls for 
Solid.genvbar, gennvv, genaicm, solve, calcpress and calcforces. Uses given freestream parameters and outputs time duration report if echo==True. Applies PG correction factor beta'''
paneller.Solid.addorder.__doc__='''addorder(self, queue, colmat, ind1, ind2): adds order for multiprocessing queue for subprocess_genaicm function to calculate AIC matrix\'s columns ind1 to ind2'''
paneller.Solid.panel_calcSn.__doc__='''panel_calcSn(p): recalculates area and normal vector for panel p'''
paneller.Solid.PG_apply.__doc__='''PG_apply(beta, a, b): applies PG correction to all panel and line positions (including wake panels) according to correction factor beta and freestream parameters a and b'''
paneller.Solid.PG_remove.__doc__='''PG_remove(beta, a, b): reverses effects from Solid.PG_apply(beta, a, b) and converts AIC matrixes and self-influence matrixes to compressible values'''

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
Returns aflpts, array of shape (2*xdisc+1, 2). Argument sweep compresses the airfoil's y axis dimension so as to apply sweep effect on thickness'''
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
utils.Mtostream.__doc__='''Mtostream(a, b): returns a matrix for coordinate system conversion to a stream-based coordinate system'''
utils.Mstreamtouni.__doc__='''Mstreamtouni(a, b): returns the inverse of utils.Mtostream(a, b)'''
utils.PG_xmult.__doc__='''PG_xmult(beta, a, b): returns a matrix to multiply point arrays so as to apply PG correction over their location, based on correction factor beta=sqrt(1.0-M**2)'''
utils.PG_inv_xmult.__doc__='''PG_inv_xmult(beta, a, b): returns the inverse of PG_xmult'''
utils.PG_vtouni.__doc__='''PG_vtouni(beta, a, b): returns matrix to multiply velocities by so as to convert them from incompressible PG circunstances to compressible values'''

control.__doc__='''
Module for definition of controls, control axes and degrees of freedom for aicraft and wing classes.

=========
functions
=========

z_rotation_matrix

=======
classes
=======

control_DOF: class to contain a degree of freedom controllable from aircraft class methods, to control deflection in one or more control axis
control_axis: class contaning a control axis and its rotation method, to rotate panels in wing class and Solid class
control: a general control definition, containing an axis, a correspondent degree of freedom and a multiplier for the given DOF to be applied to control deflection'''

control.z_rotation_matrix.__doc__='''z_rotation_matrix(th): to return a rotation matrix of angle th (radians) around z axis'''
control.control_DOF.__doc__='''
Class to model a given degree of freedom, as called from aircraft class instance

=========
variables
=========

state: numerical value for control deflection, controlled from aircraft instance

=======
methods
=======

__init__
'''
control.control_DOF.__init__.__doc__='''__init__(): initializes control_DOF instance with zero deflection'''
control.control_axis.__doc__='''
Class to define an axis around which to rotate controlled panels

=========
variables
=========

control_rot_func: lambda to compute rotation of a point around given control axis

=======
methods
=======

__init__
'''
control.control_axis.__init__.__doc__='''__init__(p0=np.array([0.0, 0.0, 0.0]), p1=np.array([0.0, 1.0, 0.0])): instantiate a control axis defined between points p0 and p1'''
control.control.__doc__='''
Class to be added to wing quadrants and controlled from aircraft instance

=========
variables
=========

DOF: degree of freedom (control_DOF instance) to be controlled externally, from aircraft instance
axis: control_axis instance
multiplier: multiplier for control deflections
paninds: index of panels rotated by given control

=======
methods
=======

__init__
addpanels
'''
control.control.__init__.__doc__='''__init__(DOF=None, p0=np.array([0.0, 0.0, 0.0]), p1=np.array([0.0, 1.0, 0.0]), multiplier=1.0): instantiate control
with controlling DOF in DOF kwarg, with axis defined between points p0 and p1, with multiplier in multiplier kwarg'''
control.control.addpanels.__doc__='''addpanels(panlist): adds panels listed in index list panlist to control.paninds list of rotated panels\'s indexes'''

wing.__doc__='''
Module containing classes for wing definition

=======
classes
=======

wing_section: section of a wing, with its airfoil points. Wing quadrants are defined between two of them
wing_quadrant: space between two wing sections
wing: a set of wing quadrants abuted to each other
'''

wing.wing_section.__doc__='''
Class containing info from a wing section

=========
variables
=========

points: array of points in section (shape 2*xdisc+1, 3)
CA_position: position of aerodynamic center (quarter-chord point) of the section at hand
c: chord of section at hand
closed: whether or not the section should be considered a closed wingtip
inverted: whether or not the airfoil in question has been inverted (as done with horizontal stabilizers)
correction: a sectional viscous correction for the section at hand (check reference on xfoil_visc module for more information)
controls: list of control objects
control_multipliers: list of multipliers for control objects
control_ptinds: list of lists of indexes of points in wing_section.points array that are included in each control\'s influence
alphas: lambda to set local angle of atack based on section\'s inviscid lift coefficient
Cls: lambda to set local viscous lift coefficient variation based on section\'s inviscid lift coefficient
Cds: lambda to set local viscous drag coefficient variation based on section\'s inviscid lift coefficient
Cms: lambda to set local viscous moment coefficient variation based on section\'s inviscid lift coefficient

=======
methods
=======

__init__
hascorrection
getcorrection
addcontrols
hascontrol
hasdeflection
getinthick
applycontrols
'''
wing.wing_section.__init__.__doc__='''__init__(afldir='', c=1.0, incidence=0.0, gamma=0.0, CA_position=np.array([0.0, 0.0, 0.0]), afl='n4412', \
        header_lines=1, xdisc=10, remove_TE_gap=True, inverse=False, xstrategy=lambda x: (np.sin(pi*x-pi/2)+1)/2, closed=False, \
            correction=None, Re=2e6, sweep=0.0):
            initializes wing section with airfoil \'afl\' collected from directory afldir (if len(afldir)==0: afldir=os.getcwd()) through utils.read_afl,
            with arguments header_lines, disc=xdisc, remove_TE_gap, inverse and strategy=xstrategy. closed kwarg sets wing section to closed wingtip if True.
            if correction!=None, a correction from module xfoil_visc is provided, and Reynolds Re should be considered to generate a viscous polar for the given section.
            Argument sweep compresses the airfoil's y axis dimension so as to apply sweep effect on thickness'''
wing.wing_section.hascorrection.__doc__='''hascorrection(): returns whether or not xfoil viscous correction is available (i. e. self.correction!=None)'''
wing.wing_section.getcorrection.__doc__='''getcorrection(Re=2e6): sets lambas self.alphas, self.Cls, self.Cds, self.Cms for wing_section instance'''
wing.wing_section.addcontrols.__doc__='''addcontrols(controls=[], control_multipliers=[], control_axpercs=[]): add controls with percentages in chordwise direction and multipliers
listed in kwarg lists, and their control objects in controls kwarg list'''
wing.wing_section.hascontrol.__doc__='''hascontrol(): return hasattr(\'controls\'), i. e. boolean on whether the section is subject to any control object'''
wing.wing_section.getinthick.__doc__='''getinthick(eta=0.5, x=0.0): get point in section in x position in x-axis and at an eta percentage of the airfoil\'s thickness (positive on
positive orientation of z axis)'''
wing.wing_section.applycontrols.__doc__='''applycontrols(ths, control_inds): applies control deflections in deflection list th (a single float value can also be provided
for deflecting all controls by the same value) by deflecting the panels listed for each control object by moving the points indexed in section by indexes present in sublists
of list of lists control_inds'''

wing.wing_quadrant.__doc__='''
Class containing information and methods for a wing quadrant, defined as the space between to wing sections

=========
variables
=========

sld: Solid instance upon whitch to compose the quadrant\'s panel patches
sect1: left (or upper, if a vertically oriented quadrant) wing section
sect2: right (or lower, if a vertically oriented quadrant) wing section
wakecomb: list of lists with panel indexes (respectively upper and lower surface panel indexes) with which to generate the wake on a given panel strip
panstrips_extra: upper wing surface\'s panel indexes, organized in sublists that contain a given chordwise wing strip
panstrips_intra: analogous list of lists for the quadrant\'s lower surface
acft: aircraft instance to which the quadrant is bound
controls: if present (check reference for wing_quadrant.hascontrol()), a dictionary containing panel objects in quadrant with keys correspondent to names set in aircraft instance
sect1_control_indlist: control panel index matrix as in class wing_section, for section 1 (left/upper)
sect2_control_indlist: control panel index matrix as in class wing_section, for section 2 (right/lower)
extraright, extraleft, intraright, intraleft_lines: lists of line unsigned indexes identifying which elements are positioned in the quadrant\'s extremities
extraright, extraleft, intraright, intraleft_points: lists of point arrays identifying which elements are positioned in the quadrant\'s extremities
ys, cs, Cls, Cms, Cds, Cls_corrected, Cds_corrected, Cms_corrected: geometrical position of measured panel strips and correspondent sectional adimensional coefficients

=======
methods
=======

__init__
set_aircraft
trim_bybody
plot_input
hascontrol
patchcompose
calc_reference
close_tip
calc_coefs
hascorrection
calc_corrected_forces
'''
wing.wing_quadrant.__init__.__doc__='''__init__(sld, sect1=None, sect2=None, control_names=[], control_axpercs_x=[], control_axpercs_thickness=[], \
        control_multipliers=[]): constructor for wing_quadrant class. sect1 and sect2 set correspondent sections in the wing quadrant (by convention:
        sect1 designates the upper/left section in the quadrant, and sect2, the lower/right section). sld binds the quadrant to the given Solid instance.
        control_names is a list of strings that designates names to controls defined in quadrant. control_axpercs_x and control_axpercs_thickness are list of lists
        that define the x and z positions in sect1 and sect2 for the axis\'s extremities: e. g. control_names=['aileron'], control_axercs_x=[[0.7, 0.75]], 
        control_axpercs_thickness=[[0.5, 0.2]] define a control called aileron whose axis should start at 70 percent chordwise in section 1 and 75 in section 2, and 
        50 percent thickness on z axis in section 1 and 20 percent in section 2. control_multipliers defines a multiplier for the axis deflection in the given quadrant,
        in respect to the axis\' DOF variable instantiated within an aircraft object. Check references for control module for further detail. Default multiplier set to 1'''
wing.wing_quadrant.set_aircraft.__doc__='''set_aircraft(acft): sets wing_quadrant.acft to aircraft instance acft in argument'''
wing.wing_quadrant.trim_bybody.__doc__='''trim_bybody(contactbody, sectside=2, tolerance=0.00005): drags points in section 1 (setside==1, upper/left) or 2 (setside==2, lower/right)
until their contact with a body (contactbody) class instance is obtained. Geometrical tolerance is set by the accordingly named kwarg. A warning is issued if any point
fails to find a point of contact with the given body, and body.patchcompose functions should be executed only after the issue is solved'''
wing.wing_quadrant.plot_input.__doc__='''plot_input(fig=None, ax=None, show=False, xlim=[], ylim=[], zlim=[], colour='blue'): plot a wireframe for the given wing quadrant, with optional axis
limits, to a previous figure and set of axes from matplotlib or a newly instantiated one (if ax==None or fig==None). if show=True, plt.show() is ran within this function. Else it should be ran externally'''
wing.wing_quadrant.hascontrol.__doc__='''hascontrol(): return whether the quadrant in question has any control (i. e. hasattr(self, 'controls'))'''
wing.wing_quadrant.patchcompose.__doc__='''patchcompose(prevlines={}, strategy=lambda x: (np.sin(pi*x-pi/2)+1)/2, ldisc=20, tolerance=5e-5): composes patches of panels for bound Solid instance, with y-axis discretization
ldisc set as np.interp(strategy(np.linspace(0.0, 1.0, ldisc)), np.array([0.0, 1.0]), np.array([sect1.CA_position[1], sect2.CA_position[1]])). Geometric tolerance for Solid.addpatch is set as kwarg tolerance'''
wing.wing_quadrant.calc_reference.__doc__='''calc_reference(axis=1): calculate wing planar surface of wing quadrant (on xy plane if axis==1 and xz if axis==2) and $\int^{sect2}_{sect1}c^2dy$, returned in tuple in given order'''
wing.wing_quadrant.close_tip.__doc__='''close_tip(sectside=2): closes wing tip at section 1 (sectside==1, upper/left) or section 2 (sectside==2, lower/right) by enforcing finite circulation around tip vortex lines. Called from
wing instance\'s patchcompose function'''
wing.wing_quadrant.calc_coefs.__doc__='''calc_coefs(alpha=0.0, axis=1): calculate sectional coefficients (self.ys, self.cs: geometric properties per section; self.Cls, self.Cds, self.Cms: actual coefficients per section)
in wing\'s axis (if axis==2, note that alpha should be interpreted as a sideslip angle). Viscous corrections, if present, are readily applied (check self.Cls_corrected, self.Cds...)'''
wing.wing_quadrant.hascorrection.__doc__='''hascorrection(): returns whether both wing quadrant sections have viscous corrections for application'''
wing.wing_quadrant.calc_corrected_forces.__doc__='''calc_corrected_forces(): calculate corrections for sectional forces based on each edge section\'s available viscous corrections.
Returns variations in coefficients return dCX, dCY, dCZ, dCl, dCm, dCn integrated along wing quadrant'''

wing.wing.__doc__='''
Class to model a whole wing through calls to methods of several queued wing quadrants

=========
variables
=========

sld: Solid instance to which the wing is bound
coefavailable: whether or not sectional coefficients have been calculated in the extent of the wing
wingquads: list of wing quadrants, from left/up to right/down
axis: the axis to which the wing is closest to allign to
acft: aircraft object to which the wing is bound
extraright, extraleft, intraright, intraleft_lines: lists of line unsigned indexes identifying which elements are positioned in the wing\'s extremities
extraright, extraleft, intraright, intraleft_points: lists of point arrays identifying which elements are positioned in the wing\'s extremities
paninds: list with all panel indexes belonging to the wing\'s patches
ys, cs, Cls, Cms, Cds, Cls_corrected, Cds_corrected, Cms_corrected: geometrical position of measured panel strips and correspondent sectional adimensional coefficients

=======
methods
=======

__init__
set_aircraft
patchcompose
trim_bybody
calc_coefs
calc_reference
close_tip
plot_input
genwakepanels
hascorrection
calc_corrected_forces
'''
wing.wing.__init__.__doc__='''__init__(sld, wingquads=[]): instantiate wing with its list of wing quadrants (listed left/up to right/down) and its bound solid (sld)'''
wing.wing.set_aircraft.__doc__='''set_aircraft(acft): sets the wing\'s bound aircraft to acft'''
wing.wing.patchcompose.__doc__='''patchcompose(ystrategy=lambda x: x, ydisc=20, tolerance=5e-5): compose patches of each wing quadrant in the wing; if ydisc is an integer, discretization is divided
along the wing according to each quadrant\'s length. Rule yspacing=np.interp(ystrategy(np.linspace(0.0, 1.0, ydisc)), np.array([0.0, 1.0]), np.array([y1, y2])) is used to set spacing of panels across the span.
If ydisc is a list, it is trimmed by reproducing the first element to the lngth of wing.wingquads, and each element is attributed to the discretization of a wing quadrant (leftmost to rightmost).
tolerance kwarg is the geometric tolerance for Solid.addpatch() method'''
wing.wing.trim_bybody.__doc__='''trim_bybody(contactbody, sectside=2, tolerance=5e-5): trim the leftmost/highest wing quadrant (sectside==1) or the rightmost/lowest wing quadrant (sectside==2) so that the wing is limited by contactbody. Check reference
for wing_quadrant.trim_bybody() for further information'''
wing.wing.calc_coefs.__doc__='''calc_coefs(alpha=0.0, axis=1): alculate sectional adimensional coefficients along the wing by calling wing_quadrant.calc_coefs() for each and concatenating resulting vectors.
Check its reference for further detail'''
wing.wing.calc_reference.__doc__='''calc_reference(): calculates S, mac, b (respectively returned in tuple) of whole wing based on call to wing_quadrant.calc_reference() for each of the component quadrants'''
wing.wing.close_tip.__doc__='''close_tip(sectside=2): impose limited circulation in wingtip vortex line segments for leftmost/highest section (sectside==1) or rightmost/lowest section (sectside==2)'''
wing.wing.plot_input.__doc__='''plot_input(fig=None, ax=None, show=False, xlim=[], ylim=[], zlim=[], colour='blue'): plots the wing as a wireframe in axes ax and figure fig. New ones are instantiated if fig==None or ax==None.
If show is True, plt.show() is ran automatically; otherwise it must be ran externally'''
wing.wing.genwakepanels.__doc__='''genwakepanels(offset=1000.0, a=0.0, b=0.0): generates wake panels for each wing quadrant, considering freestream parameters and length for wake provided in kwargs'''
wing.wing.hascorrection.__doc__='''hascorrection(): checks whether all wing quadrants in the wing have been provided with viscous corrections. It is a condition for the application of those corrections
in coefficient calculations'''
wing.wing.calc_corrected_forces.__doc__='''calc_corrected_forces(): returns dCX, dCY, dCZ, dCl, dCm, dCn considering viscous corrections'''

body.__doc__='''
Module containing definitions necessary for non-lifting body definitions

=========
functions
=========

Re2e5
Blausius_Cf_l
Prandtl_1_7th
circdefsect
squaredefsect
smooth_angle_defsect_function
standard_body
prevline_organize

NOTE: a defsect designates a function returning a body section according to arguments (y_expand, z_expand, R, center).

=======
classes
=======

body_section
body_panel
body
'''
body.Re2e5.__doc__='''Re2e5(Re): simplest possible turbulence criterion (returns Re>2e5)'''
body.Blausius_Cf_l.__doc__='''Blausius_Cf_l(Re): returns friction coefficient according to Blausius\'s laminar boundary layer formulation'''
body.Prandtl_1_7th.__doc__='''Prandtl_1_7th(Re): returns friction coefficient for laminar boundary layers according to Prandtl\'s 1-7th power rule'''
body.prevline_organize.__doc__='''prevline_organize(queue, nlines, prevlateral=[], intra=False, right=False): organizes previous lines for list to be provided to prevlines dictionary in Solid.addpatch 
method in body patch composing. Takes arguments:
* queue: a list of abuted lifting surfaces, with crescent x position of aerodynamic center of abuted section
* nlines: list with number of lines in each of the intervals between abuted lifting surfaces (e. g. [20, 10, 30] could be provided for queue=[leftwing, horizontal_empenage]). Has an orientation
opposite to \'queue\', in the x axis
* prevlateral: previous patch edge, in vortex line unsigned indexes, to fetch lines from when composing edge for new patch. If empty, new lines are created for that purpose and no previous patch is considered
in the body
* intra: to fetch lines from abuted surfaces\' lower surface (boolean)
* right: to fetch lines from queue[i].intraright_lines and queue[i].extraright_lines, if true, for each element of the provided lifting surface queue'''
body.circdefsect.__doc__='''circdefsect(y_expand=1.0, z_expand=1.0, R=1.0, center=np.array([0.0, 0.0, 0.0]), cubic=True, disc=360): returns circular section with given geometrical parameters (and
cubic polar rule, if cubic==True) for body composition'''
body.squaredefsect.__doc__='''circdefsect(y_expand=1.0, z_expand=1.0, R=1.0, center=np.array([0.0, 0.0, 0.0]), cubic=True, disc=360): returns square section with given geometrical parameters (and
cubic polar rule, if cubic==True) for body composition'''
body.smooth_angle_defsect_function.__doc__='''smooth_angle_defsect_function(r_1x=0.5, r_2x=0.5, r_1y=0.5, r_2y=0.5, ldisc=30, thdisc=20): returns a defsect (a lambda function returning 
a body section according to arguments (y_expand, z_expand, R, center)) with ellipsoidally concordant corners of semi-axis lengths r_1x, r_1y 
(lower corners\' semi-axis lengths/R, R being the section\'s largest dimension) and r_2x, r_2y (same for upper corners).
Check monoplane.py code (in LovelacePM package folder: import LovelacePM.monoplane) for an example of its use'''
body.standard_body.__doc__='''standard_body(sld, defsect=circdefsect, nose_loc=np.array([0.0, 0.0, 0.0]), nose_length=0.1, nose_thdisc=10, body_length=1.0, \
    body_width=0.1, tailcone_length=0.2, body_thdisc=60, tolerance=0.00005, nose_lift=0.0, tail_lift=0.0, z_expand=1.0, \
y_expand=1.0):
returns a body bound to solid sld of:
* straight body of length body_length, defined with ;
* input radial discretization body_thdisc for the body;
* ellipsoidal nose of length nose_length, with tip raised by nose_lift in z axis;
* tailcone of length tailcone_length, with tip raised by tail_lift in z axis;
* shape defined by defsect (a lambda function returning a body section according to arguments (y_expand, z_expand, R, center)) and scaling factor kwargs z_expand and y_expand;
* geometrical tolerance in body_panel class creation is set by kwarg tolerance.'''

body.body_panel.__doc__='''
Class with information about a panel that delimits an input body (for geometric, and not aerodynamic purposes).

=========
variables
=========

points: array (shape (3, 4)) with four points that compose the panel
center: center of the panel
Mtosys: basis conversion matrix for a panel-oriented system (z-axis normal to plane, pointing out of the solid)
Mtouni: its inverse
locpoints: points converted to local coordinate system (Mtosys@(points-center))

=======
methods
=======

__init__
'''
body.body_panel.__init__.__doc__='''__init__(p1, p2, p3, p4, tolerance=0.00005): constructor for a body panel,
recieving its four points and a geometric tolerance (tolerated areas should have at least tolerance**2)'''

body.body_section.__doc__='''
Class computing information about a section of a body in the yz plane. Should be lined, in body class instances, in the direction of the x axis.

=========
variables
=========

R: biggest dimension
coords: coordinates of points (array of shape (npoints, 3))
thetas: array with polar coordinates of points (positive as counterclockwise around x axis, on yz plane. 0 on z axis, ranging from -pi to pi)
Rs: correspondent radii, used along with array thetas
center: the sections center, around which polar coordinates will be calculated
polar_rule: lambda returning R with th=theta (polar coordinate as described before) argument. Defined as spline or linear interpolation, according to argument in constructor

=======
methods
=======
__init__
__call__
'''
body.body_section.__init__.__doc__='''__init__(center=np.array([0.0, 0.0, 0.0]), coords=np.vstack((np.sin(np.linspace(0.0, 2*pi, 360)), \
        np.cos(np.linspace(0.0, 2*pi, 360)))).T, R=1.0, y_expand=1.0, z_expand=1.0, cubic=True): constructor of a section with coordinates based on threedimensionalization
        of yz-plane coordinates coords, with maximum dimension R, around center. Scaled in y and z axis by multipliers y_expand and z_expand. If cubic is True, the object\'s
        polar_rule lambda function (responsible for interpolating points along the section, with polar coordinates as arguments) 
        will be defined as a cubic spline. Otherwise, as a linear interpolation'''
body.body_section.__call__.__doc__='''__call__(th): returns point around center, with polar coordinates in [-pi; pi], clockwise around x axis with zero on z axis, and radius as defined by
body_section.polar_rule'''

body.body.__doc__='''
Class to agregate information about a non-lifting body. Patches are composed, for it, with 4 patches angularly distributed around the body and abutted lifting surfaces between them.

=========
variables
=========

sld: Solid instance to which the body is bound
tolerance: geometric tolerance for patchcompose functions
sect_xpos: vector including x position of all section centers
sections: list of sections that compose the non-lifting body
center: \'lambda x:\' function returning the position in 3D space of the center of a transversal fuselage cut at position x in x axis
y_expand_rule, z_expand_rule: linearly interpolations defining y_expand and z_expand section arguments as a function of x
last_x: interpolation object identifying the x position of the last section (previous in positive x-axis listing)
next_x: interpolation object identifying the x position of the next section (next in positive x-axis listing)
body_panels: list with all body panel instances
paninds: patch panel index matrix for the 4 patches in a non-lifting body

=======
methods
=======

__init__
set_aircraft
find_body_intersect
plot_input
surfinterp
side_separate
line_surfinterp
theta_queueident
patchcompose
apply_eqflatplate
bodypanel_plotnormals
'''
body.body.__init__.__doc__='''__init__(sld, sections=[], tolerance=0.00005):
creates body object with sections in kwarg list, with given geometric tolerance for patch definition and bound to Solid object sld'''
body.body.set_aircraft.__doc__='''set_aircraft(acft):
binds the given non-lifting body to aircraft class instance acft'''
body.body.find_body_intersect.__doc__='''find_body_intersect(p, u, tolerance=0.00005):
find a contact point with a line defined by point p and vector u. Geometric tolerance is set for maximum distance between 
a panel and a point to define the latter\'s belonging to it'''
body.body.plot_input.__doc__='''plot_input(fig=None, ax=None, show=False, xlim=[], ylim=[], zlim=[], colour='gray'):
plot input body geometry as a wireframe. Use given figure and axes if not None, or create
new ones otheriwse. If show is True, plt.show() will be automatically ran. Otherwise, it must be externally executed'''
body.body.surfinterp.__doc__='''surfinterp(th, x):
returns point array, in 3D space, of a point in the non-lifting body\'s surface in x position in the x axis and th polar coordinate (in [-pi; pi], 0 at
z axis, positive orientation as clockwise around x axis)'''
body.body.side_separate.__doc__='''side_separate(leftqueue=[], rightqueue=[], upqueue=[], lowqueue=[], xstrategy=lambda x: (np.sin(pi*x-pi/2)+1)/2, xdisc=100):
returns x-axis spacing (given by np.interp(xstrategy(np.linspace(0.0, 1.0, xdisc+1)), np.array([0.0, 1.0]), np.array([x1, x2]))) 
and the number of lines occupying each x-axis interval between queued lifting surfaces.
returns lnlines, lxdisc, rnlines, rxdisc, unlines, uxdisc, dnlines, dxdisc. inlines designates the number of lines in intervals not occupied by lifting surfaces of queue i 
(i being d for low, u for upper, r for right, l for left queue). ixdisc indicates an array giving the x-spacing control points for queue i'''
body.body.theta_queueident.__doc__='''theta_queueident(queue, xspacing, intra=False, right=False, queueident=\'l\'): 
sets polar coordinates (in [-pi; pi], 0 at z axis, positive orientation as clockwise around x axis) for points defined by xspacing array, along a queue, 
following upper surface points (if intra==False, lower if otherwise) of surfaces in queue list. 
Points in the rightmost sections of the wings contained in queue are fetched if right is True, or from the leftmost section if otherwise.
queueident identifies the requested queue (d for low, u for upper, r for right, l for left queue)'''
body.body.patchcompose.__doc__='''patchcompose(leftqueue=[], rightqueue=[], upqueue=[], leftqueue=[], xstrategy=lambda x: x, xdisc=100, thstrategy=lambda x: x, thdisc_upleft=20, thdisc_upright=20, \
    thdisc_downleft=20, thdisc_downright=20, tolerance=5e-5): compose four panel patches (upleft: between upper and left queues; downright: between lower and rightmost queues; etc.) to model the body at hand.
    thdisc_[patch] arguments set discretization levels for each of the four patches. xdisc sets discretization in the x axis, as np.interp(xstrategy(np.linspace(0.0, 1.0, xdisc+1)), np.array([0.0, 1.0]), 
    np.array([x1, x2]))). Angular discretization in patches is defined likewise. tolerance is provided as an argument for Solid.addpatch() routine (check its referrence for further detail)'''
body.body.apply_eqflatplate.__doc__='''apply_eqflatplate(rho=1.225, Uinf=1.0, mu=1.72*10e-5, turbulent_criterion=Re2e5, Cf_l_rule=Blausius_Cf_l, Cf_t_rule=Prandtl_1_7th): applies, 
according to given turbulence criteria (function of local Reynolds returning boolean, True if turbulent) and rules for turbulent and laminar friction coefficients (function of local Reynolds,
returning floats), values for sld.Cfs[p] in the Solid instance the body at hand is bound to. Leads to approximation of viscous effects on the non-lifting body. Custom rules may be applied by means of kwargs'''
body.body.bodypanel_plotnormals.__doc__='''bodypanel_plotnormals(xlim=[], ylim=[], zlim=[], factor=1.0): function to plot normal vectors (given by third row of body_panel.Mtosys coordinate system matrix)
of each of the body panels defining the body at hand. Used for debugging input geometries. Works as Solid.plotnormals, the referrence of which should be checked for further detail'''

aircraft.__doc__='''
Module with classes to process an aircraft in its entirety, with fully automated case studies

=======
classes
=======

aircraft
'''
aircraft.aircraft.__doc__='''
Class to contain information on a whole aircraft and its case study, with an attached Solid instance and the due input wings and bodies

=========
variables
=========

wings: list of lifting surfaces in case study
bodies: list of non-lifting surfaces in case study
CG: center of gravity of the aircraft at hand
CX, CY, CZ, CL, CD, Cm, Cn, Cl: adimensional force coefficients for the whole configuration, with moments centered about the given CG
dCX, dCY, ...: variation in adimensional force coefficients due to viscous corrections
Sref, cref, bref, AR: reference quantities for calculation of adimensional coefficients
a, b, p, q, r, Uinf: freestream parameters
plotlim: maximum dimension of the aircraft stretching from the origin of its coordinate system, used as default for plot_input() and plotgeometry() limits
controlset: dictionary reporting all control degrees of freedom for the aircraft, as generated by wing quadrant constructors
forcesavailable, stabavailable, massavailable, hascorrections: whether, respectively, forces, stability derivatives, mass information and viscous corrections are available for the aircraft at hand
stabderivative_dict: a dictionary of dictionaries containing last calculated inviscid stability derivatives. E. G. Key \'a\' given to the first dictionary layer points to a dictionary with CXa, CYa...
Key \'dCl\', for example, given to the \'a\' subdictionary, points to Cla[alpha] stability derivative (stabderivative_dict[\'a\'][\'dCl\'])

=======
methods
=======

__init__
parameter_report
balance
transp_byvec
transp_to_cg
hascontrol
addwake
calcforces
freestream_derivatives
calcstab
stabreport
edit_parameters
bodies_eqflatplate_apply
eulersolve
forces_report
plot_input
plotgeometry
'''
aircraft.aircraft.__init__.__doc__='''__init__(sld, elems=[], Sref=0.0, cref=0.0, bref=0.0, echo=True, CG=np.array([0.0, 0.0, 0.0])): aircraft constructor, using bound Solid object as input, and all
wings and bodies in elems list. Reference quantities, if provided as zero/not provided, will be calculated from the first input wing if any is present (or set to 1.0 if none is provided). If echo
is True, aircraft.parameter_report() will be ran. The function will run silently, otherwise. Variable aircraft.CG will be set in place. All Xavailable variables will initially be set to False'''
aircraft.aircraft.parameter_report.__doc__='''parameter_report(): reports in prompt the freestream parameters and control deflections (edited by aircraft.edit_parameters()) set for the case study at hand.
Angles are input and output in degrees and kept in variables as radians. Angular velocities are normalized by current Uinf and reference dimensions before being input and output'''
aircraft.aircraft.balance.__doc__='''balance(SM=0.1, echo=True): transports the aircraft\'s center of gravity based on previously calculated stability derivatives, or calculates them in place
if none is available. If echo is True, previous and newly set static margins are printed on stdio, and occasionally necessary calls to method aircraft.calcstab are ran with kwarg echo=True (non-silently)'''
aircraft.aircraft.transp_byvec.__doc__='''transp_byvec(vec): shifts the aircraft\'s center of gravity by a vector \'vec\' in 3D space. Recalculates adimensional moment coefficients, stability derivatives
and viscosity corrections to moment coefficients according to the shift performed'''
aircraft.aircraft.transp_to_cg.__doc__='''transp_to_CG(CG): transports the center of gravity of the aircraft at hand to point in kwarg CG, through successive calls to aircraft.transp_byvec()'''
aircraft.aircraft.hascontrol.__doc__='''hascontrol(): evaluates whether or not the aircraft has any control surface, based on the length of its controlset dictionary. The dictionary is, in turn, generated 
from each wing quadrant\'s control lists'''
aircraft.aircraft.addwake.__doc__='''addwake(offset=1000.0, wakedisc=1, strategy=lambda x: ((np.exp(x)-1.0)/(exp(1)-1.0))**2): 
generates wakes for all lifting surfaces based on a given offset from the trailing edge (i. e. length of the trailing edge vortexes) and currently
set freestream parameters. Must be ran after parameter definition and patch composing, though before Euler solution and post-processing.
the wake is discretized in wakedisc panels queued in x-axis, with offset from trailing edge set as offset*strategy(np.linspace(0.0, 1.0, wakedisc+1))'''
aircraft.aircraft.calcforces.__doc__='''calcforces(echo=True): calculates forces acting on the aircraft based on the Solid instance aircraft.sld\'s last ran Euler solution. If echo is True, method
aircraft.forces_report() is ran in place. Viscous corrections are also computed in place (see variables dCX, dCY,... in the class\'s documentation)'''
aircraft.aircraft.freestream_derivative.__doc__='''freestream_derivative(par=\'a\'): calculates derivative of unitary streamwise (u) vector and its perpendicular vector v, parallel to lift, in respect to
freestream parameters \'alpha/a\' and \'beta/b\'. Returns (du, dv), both arrays of shape (3). Called within aircraft.calcstab() method for adjoint calculation of stability derivatives'''
aircraft.aircraft.calcstab.__doc__='''calcstab(echo=True): calculates all stability derivatives for stabderivative_dict (check reference for aircraft class\'s variables) through adjoint method. If echo is 
provided as True, method aircraft.stabreport() is called in place. Otherwise, calculations are ran silently'''
aircraft.aircraft.stabreport.__doc__='''stabreport(): reports current stability derivatives, or calculates them in place if none is available'''
aircraft.aircraft.edit_parameters.__doc__='''edit_parameters(pardict, echo=True): Changes parameters in pardict to their provided values. E. g. pardict={\'a\':10.0, \'b\':5.0} sets angle of attack to 10 deg
and sideslip angle to 5 deg. If echo is True, method aircraft.parameter_report() is ran at the end. If the key to a control surface as registered in aircraft.controlset is provided, its state variable within
the controlset dictionary will be set to the value provided in pardict, in degrees. Angles are input and output in degrees and kept in variables in radians. Angular velocities are provided normalized by Uinf 
and reference dimensions'''
aircraft.aircraft.bodies_eqflatplate_apply.__doc__='''bodies_eqflatplate_apply(rho=1.225, mu=1.72*10e-5, turbulent_criterion=Re2e5, Cf_l_rule=Blausius_Cf_l, Cf_t_rule=Prandtl_1_7th): applies equivalent flat
plate correction for friction coefficient according to given criteria for all bodies bound to the aircraft according to current freestream parameters. Check body.apply_eqflatplate() reference for further detail'''
aircraft.aircraft.eulersolve.__doc__='''eulersolve(echo=True, damper=0.0, wakeiter=0): runs Solid.eulersolve() with the present freestream parameters set for the aircraft. Time report is presented if echo is True.
A Thikhonov regularization (with damping coefficient set by damper kwarg) can be performed with the AIC matrix if damper is non-zero. wakeiter iterations of wake rollup are executed'''
aircraft.aircraft.forces_report.__doc__='''forces_report(): prints all force and moment coefficients available (and their viscous corrected equivalents). If none is available, method aircraft.calcforces() is ran'''
aircraft.aircraft.plot_input.__doc__='''plot_input(xlim=[], ylim=[], zlim=[]): plots all elements of the aircraft as wireframes for input geometry debugging. Sets aircraft.plotlim as plot limits if none
is provided'''
aircraft.aircraft.plotgeometry.__doc__='''plotgeometry(xlim=[], ylim=[], zlim=[], factor=1.0, velfield=True): plots the aircraft\'s panelling as a wireframe. If velfield is set to True, local velocities are
also drawn as vectors (multiplied by scaling factor kwarg for easier visualization)'''

aerodynamic_output.__doc__='''
Module containing functions for plotting aerodynamic results

=========
functions
=========

plot_Cps
plot_Cls
plot_gammas
plot_Cds
plot_Cms
'''
aerodynamic_output.plot_Cps.__doc__='''plot_Cps(sld, elems=[], xlim=[], ylim=[], zlim=[-2.5, 2.5]): plots Cps in panels belonging to elements (wings and bodies) in elems list according to Cps available in
bound Solid instance sld, using a 3D scattergram. 3D plotting limits are set using kwarg lists. Default z-axis (Cp) limits for plotting are set by default to interval [-2.5; 2.5]'''
aerodynamic_output.plot_Cls.__doc__='''plot_Cls(sld, alpha=0.0, wings=[], axis=1): plots sectional lift coefficients accross wings in kwarg list, based on given angle of attack and the axis along which 
the lifting surfaces stretch (1 if y axis, wings, 2 if z axis, vertical fins. alpha should be interpreted as sideslip angle if axis==2). Bound Solid instance from which to gather Cp info should be provided
in arg sld. Viscous and inviscid calculations for coefficients will be contrasted'''
aerodynamic_output.plot_gammas.__doc__='''plot_gammas(sld, alpha=0.0, Uinf=1.0, wings=[], axis=1): plots circulation along wing sections using info as provided for function aerodynamic_output.plot_Cls(),
with a given freestream veloctity. Check help(plot_Cls) for further information on argument structure'''
aerodynamic_output.plot_Cds.__doc__='''plot_Cds(sld, alpha=0.0, wings=[], axis=1): plots sectional drag coefficients accross wings in kwarg list, based on given angle of attack and the axis along which 
the lifting surfaces stretch (1 if y axis, wings, 2 if z axis, vertical fins. alpha should be interpreted as sideslip angle if axis==2). Bound Solid instance from which to gather Cp info should be provided
in arg sld. Viscous and inviscid calculations for coefficients will be contrasted'''
aerodynamic_output.plot_Cms.__doc__='''plot_Cms(sld, alpha=0.0, wings=[], axis=1): plots sectional moment coefficients accross wings in kwarg list, based on given angle of attack and the axis along which 
the lifting surfaces stretch (1 if y axis, wings, 2 if z axis, vertical fins. alpha should be interpreted as sideslip angle if axis==2). Bound Solid instance from which to gather Cp info should be provided
in arg sld. Viscous and inviscid calculations for coefficients will be contrasted'''

xfoil_visc.__doc__='''
Module containing automation of Xfoil (by Mark Drella, MIT License) for viscous correction calculations, according to strip theory.
WARNING: depends on Xfoil being registered to PATH environment variable in the machine at hand

=======
classes
=======

polar_correction

=========
functions
=========

polar_data
read_polar

=======
classes
=======

polar_correction
'''
xfoil_visc.polar_data.__doc__='''polar_data(name='n4412', afldir='', ext_append=True, aseq=[-5.0, 20.0, 1.0], visc=True, Re=3e6, M=0.03, iter=300, flap=None, npan=300, LE_con=0.4, inverse=False): 
returns alphas, Cls, Cds, Cms arrays for Xfoil-calculated polar. Airfoil is searched for (with .dat extension if ext_apped==True) in afldir directory, or current directory if it is of length 0. 
Elements of aseq list should be, respectively, first angle of attack, last angle of attack and desired step in AOAs for polar composition. 
visc==True sets Reynolds number to kwarg Re. M sets Mach number. iter kwarg sets maximum number of iterations for quasi-simultaneous viscid inviscid coupling. 
flap, if present, defines a length=3 list with, respectively, the x-axis chord percentage for the flap\'s axis, the thickness percentage for the
given axis and its deflection in degrees. npan and LE_con are, respectively, the number of panels to be used for Xfoil\'s panel method and their concentration in the leading edge.
Chosen angles of attack are inverted (for horizontal tails and likewise inverted airfoil surfaces) if inverse is True.

WARNING: aseq and iter kwargs may have to be changed to obtain better polars. Non-converged AOAs are not returned: check alphas array in returned values to detect converged angles of attack'''
xfoil_visc.read_polar.__doc__='''read_polar(poldir='', polname='n4412', ext_append=True, echo=True, cubic=True): reads polar from directory poldir (or current directory if none is provided)
with name polname. If ext_append is True, extension.plr is expected. Return is an instance of class polar_correction. If cubic is True, polar interpolation functions will be created as cubic splines.
Otherwise, linear interpolations shall be used. If kwarg echo is True, polar data shall be output as it is read'''
xfoil_visc.polar_correction.__doc__='''
Class containing information on a given viscous polar obtained from Xfoil.

WARNING: not all polars obtained from Xfoil can be readily used, for unconvergence of its quasi-simultaneous viscous-inviscid coupling method may compromise some of the returned AOAs, which are
made unavailable for use by this class. It is advisable that viscous polars are created separately from other aircraft analysis scripts, saved using polar_correction.dump() method and loaded for
aircraft analysis using function xfoil_visc.read_polar()

WARNING: polars are linearly interpolated between available calculated Reynolds numbers so as to fit other Reynolds intervals. Re_high and Re_low variables must be set so as to achieve maximum precision
with this interpolation

.plr file structure:
[first line with number of inviscid polar lines]
Inviscid polar with alpha, Cl, Cd and Cm information, respectively, in four columns for all calculated AOAs
[arbitrary low Reynolds number value for polar, separated by simple spacing from number of collected AOAs for the chosen low Reynolds number]
Low Reynolds viscous polar with alpha, Cl, Cd and Cm information, respectively, in four columns for all calculated AOAs
[arbitrary high Reynolds number value for polar, separated by simple spacing from number of collected AOAs for the chosen high Reynolds number]
High Reynolds viscous polar with alpha, Cl, Cd and Cm information, respectively, in four columns for all calculated AOAs
EOF

=========
variables
=========
Re_low: lower Reynolds number provided for polar interpolation between Reynolds numbers
Re_high: higher Reynolds number provided for polar interpolation between Reynolds numbers
alphas_Re_low, Cls_Re_low, Cds_Re_low, Cms_Re_low: polar results for low Reynolds
alphas_Re_high, Cls_Re_high, Cds_Re_high, Cms_Re_high: polar results for high Reynolds
alphas_inviscid, Cls_inviscid, Cds_inviscid, Cms_inviscid: polar results for inviscid calculations
alphas_Re_low_fun, Cls_Re_low_fun, Cds_Re_low_fun, Cms_Re_low_fun: polar interpolations based on inviscid Cl input for low Reynolds
alphas_Re_high_fun, Cls_Re_high_fun, Cds_Re_high_fun, Cms_Re_high_fun: polar interpolations based on inviscid Cl input for high Reynolds

=======
methods
=======

__init__
dump
create_functions
__call__
'''
xfoil_visc.polar_correction.__init__.__doc__='''__init__(name='n4412', ext_append=True, aseq=[-10.0, 20.0, 2.0], Re_low=2e6, Re_high=3e6, Mach=0.03, flap=None, iter=300, cubic=True):
create polar correction for named airfoil, appending extension .dat for its file name if ext_append is True, between AOAs aseq[0] and aseq[1] with AOA step aseq[2]. Linear interpolation
for Reynolds number extrapolation (check reference for xfoil_visc.polar_correction class for more info) is done between Reynolds numbers Re_low and Re_high. Mach number is set to kwarg
Mach. Maximum iteration for viscous-inviscid coupling is set to kwarg iter. If cubic is True, polars will extrapolate data for AOAs using a cubic spline. A linear interpolation will be used
otherwise. flap, if present, defines a length=3 list with, respectively, the x-axis chord percentage for the flap\'s axis, the thickness percentage for the
given axis and its deflection in degrees.'''
xfoil_visc.polar_correction.dump.__doc__='''dump(poldir='', polname='n4412', ext_append=True, echo=True): prints polar in .plr format (appending extension if ext_append is True) to file polname in directory
poldir (or pwd if len(poldir)==0). Data is printed to stdio as dumped into file if echo is True'''
xfoil_visc.polar_correction.create_functions.__doc__='''create_functions(cubic=True): creates cubic (or linear if cubic==False) interpolations for variables:
alphas_Re_low_fun, Cls_Re_low_fun, Cds_Re_low_fun, Cms_Re_low_fun: polar interpolations based on inviscid Cl input for low Reynolds
alphas_Re_high_fun, Cls_Re_high_fun, Cds_Re_high_fun, Cms_Re_high_fun: polar interpolations based on inviscid Cl input for high Reynolds'''
xfoil_visc.polar_correction.__call__.__doc__='''__call__(Re=2e6, inverse=False): returns lambdas linearly interpolating low and high Reynolds polars for given kwarg Reynolds number Re, in tuple
(alphas, Cls, Cds, Cms). inverse kwarg should be set to True if an inverted airfoil (e. g. for horizontal tails) is to be used'''