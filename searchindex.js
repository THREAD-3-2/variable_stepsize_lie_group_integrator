Search.setIndex({"docnames": ["RKMK_var_step", "index", "matlab", "tredpend"], "filenames": ["RKMK_var_step.rst", "index.rst", "matlab.rst", "tredpend.rst"], "titles": ["Runge-Kutta-Munthe-Kaas (RKMK) methods with variable step size", "Documentation of <cite>variable_stepsize_lie_group_integrator</cite>", "MATLAB code", "The N-fold 3D pendulum"], "terms": {"The": [0, 1], "underli": 0, "idea": 0, "i": [0, 1, 2, 3], "express": [0, 3], "vector": [0, 2, 3], "field": [0, 2, 3], "f": [0, 2, 3], "mathfrak": [0, 3], "x": [0, 2, 3], "mathcal": [0, 3], "m": [0, 2, 3], "vert_m": [0, 3], "psi_": [0, 3], "where": [0, 3], "infinitesim": [0, 3], "gener": [0, 3], "psi": 0, "transit": [0, 3], "action": [0, 2, 3], "rightarrow": [0, 3], "g": [0, 3], "thi": [0, 2, 3], "allow": 0, "u": [0, 2, 3], "transform": 0, "problem": [0, 3], "from": [0, 1, 2, 3], "manifold": 0, "lie": [0, 2, 3], "algebra": [0, 2], "which": [0, 2], "we": [0, 2, 3], "can": [0, 3], "perform": [0, 1], "time": [0, 2, 3], "integr": [0, 1], "map": [0, 2, 3], "result": 0, "back": 0, "repeat": 0, "up": 0, "final": 0, "more": 0, "explicitli": 0, "let": [0, 3], "h_n": [0, 3], "n": [0, 1, 2], "th": [0, 2, 3], "updat": 0, "y_n": 0, "y_": 0, "1": [0, 2, 3], "begin": [0, 3], "align": [0, 3], "case": 0, "sigma": [0, 2], "0": 0, "dot": [0, 2, 3], "t": [0, 2, 3], "textrm": 0, "dexp": 0, "_": [0, 3], "circ": 0, "exp": 0, "t_": [0, 2, 3], "sigma_1": 0, "end": [0, 3], "approx": 0, "comput": [0, 1, 2, 3], "One": 0, "approach": 0, "vari": 0, "base": [0, 2], "embed": 0, "pair": [0, 1, 2], "space": [0, 3], "consist": 0, "princip": 0, "order": [0, 1, 2], "p": [0, 2], "us": [0, 2], "propag": 0, "numer": [0, 2], "solut": [0, 1, 2], "togeth": 0, "some": 0, "auxiliari": 0, "tild": 0, "onli": 0, "obtain": [0, 1], "an": 0, "estim": 0, "local": [0, 2], "error": [0, 2], "turn": 0, "deriv": [0, 2], "adjust": 0, "formula": 0, "attempt": 0, "keep": 0, "approxim": [0, 1], "equal": 0, "user": 0, "defin": [0, 2, 3], "toler": [0, 1, 2], "tol": [0, 2], "everi": 0, "both": 0, "ar": [0, 1, 3], "appli": 0, "solv": 0, "od": [0, 2, 3], "yield": 0, "two": 0, "_1": 0, "respect": 0, "same": 0, "now": [0, 2, 3], "distanc": 0, "measur": [0, 1], "between": 0, "provid": 0, "e_": 0, "truncat": 0, "thu": [0, 3], "c": 0, "h_": 0, "o": 0, "h": [0, 2], "2": [0, 2, 3], "aim": 0, "one": 0, "mai": 0, "type": 0, "theta": 0, "left": [0, 3], "frac": 0, "right": [0, 2, 3], "tfrac": 0, "typic": 0, "chosen": 0, "8": 0, "9": 0, "If": 0, "e_n": 0, "reject": 0, "henc": [0, 3], "redo": 0, "In": [0, 3], "our": [0, 3], "code": [0, 1], "experi": [0, 2], "fold": [0, 1, 2], "3d": [0, 1], "pendulum": [0, 1, 2], "compar": [0, 1], "constant": [0, 1, 3], "consid": [0, 3], "come": [0, 1, 3], "dormand": [0, 1, 2], "princ": [0, 1, 2], "dopri": [0, 1], "5": [0, 1, 2], "4": [0, 1, 2], "denot": [0, 1, 3], "set": [0, 2], "10": 0, "6": [0, 2], "system": [0, 2, 3], "scheme": [0, 2], "fix": [0, 3], "number": [0, 2], "requir": 0, "rkmk5": [0, 1, 2], "comparison": [0, 1], "occur": 0, "3": [0, 2, 3], "euclidean": 0, "norm": 0, "ambient": [0, 3], "mathbb": [0, 3], "r": [0, 2, 3], "6n": [0, 3], "connect": [0, 2, 3], "pendula": [0, 2], "qualiti": [0, 1], "against": [0, 1], "refer": [0, 1], "ode45": [0, 1, 2], "matlab": [0, 1], "strict": [0, 1], "variabl": [1, 3], "step": [1, 2], "size": [1, 2], "method": [1, 2], "rkmk": [1, 2], "also": 1, "test": 1, "exampl": 1, "part": [1, 2, 3], "sourc": [1, 2], "develop": 1, "depart": 1, "mathemat": [1, 2, 3], "scienc": 1, "ntnu": 1, "paper": 1, "celledoni": [1, 2, 3], "\u00e7okaj": [1, 2, 3], "leon": [1, 2, 3], "murari": [1, 2, 3], "owren": [1, 2, 3], "2021": [1, 3], "intern": [1, 2, 3], "journal": [1, 2, 3], "arxiv": [1, 2], "rung": [1, 2], "kutta": [1, 2], "munth": [1, 2], "kaa": [1, 2], "equat": [1, 2], "motion": 1, "src": 1, "lie_group_funct": 1, "equations_of_mot": 1, "helpful_funct": 1, "index": [1, 2], "search": 1, "page": 1, "document": 2, "sphinx": 2, "matlabdomain": 2, "extens": 2, "modul": 2, "contain": 2, "follow": [2, 3], "file": 2, "main": 2, "function": [2, 3], "adapt": 2, "ref": 2, "e": 2, "A": 2, "d": 2, "b": 2, "dynam": [2, 3], "framework": 2, "group": [2, 3], "mechan": [2, 3], "99": 2, "58": 2, "88": 2, "vecfield": 2, "paramet": 2, "hand": 2, "side": 2, "t_n": 2, "return": 2, "rkmk45": 2, "stepsiz": 2, "variablestepcomparison": 2, "z0": 2, "initi": 2, "valu": 2, "instant": 2, "new": 2, "exprodrigu": 2, "exponenti": 2, "so": 2, "input": 2, "element": 2, "repres": 2, "r3": 2, "expse3": 2, "se": [2, 3], "compon": 2, "v": 2, "correspond": 2, "skew": 2, "symmetr": [2, 3], "matrix": [2, 3], "hat": [2, 3], "translat": 2, "3x4": 2, "3x3": 2, "rotat": 2, "exponentialse3n": 2, "actionse3": 2, "actionse3n": 2, "dexpinvse3": 2, "invers": 2, "6x1": 2, "dexp_sigma": 2, "dexpinvse3n": 2, "fmanitoalgebra": 2, "q": [2, 3], "w": 2, "l": 2, "rh": 2, "posit": 2, "q1": 2, "qn": 2, "": [2, 3], "angular": [2, 3], "veloc": [2, 3], "w1": 2, "wn": 2, "qi": 2, "length": 2, "mass": [2, 3], "assemblef": 2, "precis": 2, "need": 2, "becom": 2, "here": [2, 3], "assembl": 2, "form": 2, "assemblem": 2, "inertia": 2, "block": [2, 3], "funcq": 2, "z": 2, "q2": 2, "w2": 2, "_i": [2, 3], "w_i": 2, "q_i": [2, 3], "funcw": 2, "qp": 2, "wp": 2, "initializestat": 2, "extractq": 2, "extract": 2, "extractw": 2, "It": 2, "associ": 2, "getvec": 2, "vec": 2, "term": [2, 3], "getblock": 2, "mat": 2, "j": [2, 3], "reorder": 2, "format": 2, "describ": 3, "specif": 3, "chain": 3, "whose": 3, "evolv": 3, "act": 3, "phase": 3, "present": 3, "subject": 3, "graviti": 3, "model": 3, "rigid": 3, "massless": 3, "link": 3, "serial": 3, "spheric": 3, "joint": 3, "first": 3, "point": 3, "place": 3, "origin": 3, "neglect": 3, "friction": 3, "interact": 3, "among": 3, "lee": 3, "leok": 3, "mcclamroch": 3, "2018": 3, "omit": 3, "detail": 3, "configur": 3, "m_i": 3, "euler": 3, "lagrang": 3, "q_1": 3, "q_n": 3, "subset": 3, "3n": 3, "omega_1": 3, "omega_n": 3, "kinemat": 3, "omega_i": 3, "quad": 3, "written": 3, "omega": 3, "sum_": 3, "substack": 3, "neq": 3, "m_": 3, "ij": 3, "omega_j": 3, "q_j": 3, "big": 3, "m_j": 3, "gl_i": 3, "e_3": 3, "bmatrix": 3, "r_1": 3, "vdot": 3, "r_n": 3, "ii": 3, "nm_j": 3, "l_i": 3, "2i_3": 3, "k": 3, "m_k": 3, "l_il_j": 3, "_j": 3, "ji": 3, "text": 3, "max": 3, "i_3": 3, "find": 3, "foral": 3, "sinc": 3, "linear": 3, "invert": 3, "see": 3, "a_": 3, "a_q": 3, "rewrit": 3, "h_1": 3, "a_1": 3, "a_n": 3, "r_i": 3, "a_i": 3, "h_i": 3, "given": 3, "simeq": 3}, "objects": {"": [[2, 0, 0, "-", "src"], [2, 0, 1, "", "src"]], "src.Lie_group_functions": [[2, 1, 1, "", "actionSE3"], [2, 1, 1, "", "actionSE3N"], [2, 1, 1, "", "dexpinvSE3"], [2, 1, 1, "", "dexpinvSE3N"], [2, 1, 1, "", "expRodrigues"], [2, 1, 1, "", "expSE3"], [2, 1, 1, "", "exponentialSE3N"]], "src.equations_of_motion": [[2, 1, 1, "", "FuncQ"], [2, 1, 1, "", "FuncW"], [2, 1, 1, "", "assembleF"], [2, 1, 1, "", "assembleM"], [2, 1, 1, "", "assembleR"], [2, 1, 1, "", "fManiToAlgebra"], [2, 1, 1, "", "initializeStat"]], "src.helpful_functions": [[2, 1, 1, "", "extractq"], [2, 1, 1, "", "extractw"], [2, 1, 1, "", "getBlock"], [2, 1, 1, "", "getVec"], [2, 1, 1, "", "hat"], [2, 1, 1, "", "reorder"]], "src.integrators": [[2, 1, 1, "", "RKMK45"], [2, 1, 1, "", "RKMK5"], [2, 1, 1, "", "variableStepComparison"]], "src": [[2, 2, 1, "", "main"]]}, "objtypes": {"0": "mat:module", "1": "mat:function", "2": "mat:script"}, "objnames": {"0": ["mat", "module", "MATLAB module"], "1": ["mat", "function", "MATLAB function"], "2": ["mat", "script", "MATLAB script"]}, "titleterms": {"rung": 0, "kutta": 0, "munth": 0, "kaa": 0, "rkmk": 0, "method": 0, "variabl": 0, "step": 0, "size": 0, "document": 1, "variable_stepsize_lie_group_integr": 1, "content": 1, "indic": 1, "tabl": 1, "matlab": 2, "code": 2, "src": 2, "integr": 2, "lie_group_funct": 2, "equations_of_mot": 2, "helpful_funct": 2, "The": 3, "n": 3, "fold": 3, "3d": 3, "pendulum": 3, "equat": 3, "motion": 3}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx": 56}})