import paneller
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
paneller.Solid.genwakepanels.__doc__='''genwakepanels(self, wakecombs=[], wakeinds=[], offset=1000.0, a=0.0, b=0.0): add wake panels with combinations wakecombs (list of lists, including -1s to disregard upper
or lower surface panels) with indup and indown reported in list of lists wakeinds, to be trimmed to the same length of wakecombs with default values indup=0, indown=2'''
paneller.Solid.panel_getcoords.__doc__='''panel_getcoords(self, p): return coordinates for four (or three, if triangular) extremities of Panel object instance \'p\''''
paneller.Solid.line_getcoords.__doc__='''line_getcoords(self, ind): return coordinate matrix (shape (3, 2)) for line with signed index ind. Coordinates are returned inverted if index ind
is presented as negative'''
paneller.Solid.line_getvec.__doc__='''line_getvec(self, ind): return (oriented based on signalled index ind) the vector correspondent to the vortex line segment of signalled index ind'''
paneller.Solid.line_midpoint.__doc__='''line_midpoint(self, ind): returns midpoint for line of signed index ind'''
paneller.Solid.nvect_diradjust.__doc__='''nvect_diradjust(self, patchinds, vect): adjust all nvector variables in panels whose indexes are contained in index matrix \'patchinds\' to point in a direction 
given by vector vect\'s positive orientation'''
paneller.Solid.nvect_radadjust.__doc__='''nvect_radadjust(self, patchinds, center, inwards=False): adjust all nvector variables in panels whose indexes are contained in index matrix \'patchinds\' to point
outwards (or inwards if inwards==True) from center in \'center\' argument'''