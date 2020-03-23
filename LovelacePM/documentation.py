import paneller
import utils
import control
import wing
import body
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
            correction=None, Re=2e6):
            initializes wing section with airfoil \'afl\' collected from directory afldir (if len(afldir)==0: afldir=os.getcwd()) through utils.read_afl,
            with arguments header_lines, disc=xdisc, remove_TE_gap, inverse and strategy=xstrategy. closed kwarg sets wing section to closed wingtip if True.
            if correction!=None, a correction from module xfoil_visc is provided, and Reynolds Re should be considered to generate a viscous polar for the given section'''
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

Re2e6
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
body.Re2e6.__doc__='''Re2e6(Re): simplest possible turbulence criterion (returns Re>2e6)'''
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