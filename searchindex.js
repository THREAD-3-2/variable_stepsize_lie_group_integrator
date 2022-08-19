Search.setIndex({"docnames": ["index", "installation", "matlab", "tredpend", "vi"], "filenames": ["index.rst", "installation.rst", "matlab.rst", "tredpend.rst", "vi.rst"], "titles": ["Documentation of <cite>variable_stepsize_lie_group_integrator</cite>", "Installation", "MATLAB code", "The N-fold 3D pendulum", "Variational integrator"], "terms": {"i": [0, 2, 3, 4], "matlab": 0, "code": 0, "comparison": 0, "perform": 0, "constant": [0, 3], "variabl": [0, 3, 4], "step": [0, 4], "size": 0, "method": 0, "The": [0, 4], "rkmk": 0, "pair": 0, "come": [0, 3], "from": [0, 3, 4], "dormand": 0, "princ": 0, "dopri": 0, "5": 0, "4": 0, "also": 0, "denot": [0, 3, 4], "compar": 0, "order": 0, "rkmk5": [0, 2], "ar": [0, 3, 4], "test": 0, "n": [0, 4], "fold": 0, "3d": 0, "pendulum": 0, "exampl": [0, 1], "qualiti": 0, "approxim": [0, 4], "measur": 0, "against": 0, "refer": [0, 4], "solut": 0, "obtain": [0, 4], "ode45": 0, "strict": 0, "toler": 0, "part": [0, 3], "sourc": [0, 2], "develop": [0, 1], "depart": 0, "mathemat": [0, 3], "scienc": 0, "ntnu": 0, "paper": 0, "elena": [0, 3], "celledoni": [0, 3], "ergi": [0, 3], "\u00e7okaj": [0, 3], "andrea": [0, 3], "leon": [0, 3], "david": [0, 3], "murari": [0, 3], "brynjulf": [0, 3], "owren": [0, 3], "2021": [0, 3], "intern": [0, 3, 4], "journal": [0, 3], "comput": [0, 1, 3, 4], "arxiv": 0, "instal": 0, "virtual": 0, "environ": 0, "equat": 0, "motion": 0, "variat": 0, "integr": 0, "continu": 0, "set": 0, "discret": 0, "euler": [0, 3], "lagrang": [0, 3], "src": 0, "index": 0, "modul": [0, 2], "search": 0, "page": 0, "If": 1, "you": 1, "have": [1, 4], "support": 1, "python": 1, "your": 1, "can": [1, 3, 4], "packag": 1, "like": 1, "so": 1, "creat": 1, "call": 1, "venv": 1, "python3": 1, "m": [1, 3, 4], "activ": 1, "bin": 1, "project": [1, 4], "edit": 1, "state": 1, "pip": 1, "e": [1, 4], "option": 1, "jupyterlab": 1, "jupyt": 1, "widget": 1, "ipywidget": 1, "For": [1, 4], "list": 1, "requir": [1, 4], "txt": 1, "need": 1, "r": [1, 3, 4], "To": [1, 4], "build": 1, "document": [1, 2], "local": 1, "run": 1, "cd": 1, "doc": 1, "make": 1, "doctest": 1, "check": 1, "work": [1, 4], "html": 1, "thi": [2, 3, 4], "sphinx": 2, "us": [2, 4], "matlabdomain": 2, "extens": 2, "see": [2, 3], "here": [2, 3], "se": [2, 3], "3": [2, 3, 4], "relat": 2, "function": [2, 3, 4], "contain": [2, 4], "follow": [2, 3], "file": 2, "main_const_vs_var_stepsize_comparison": 2, "rkmk45": 2, "we": 3, "describ": 3, "specif": 3, "problem": 3, "chain": 3, "connect": 3, "whose": 3, "dynam": [3, 4], "evolv": 3, "t": [3, 4], "2": [3, 4], "mechan": 3, "system": [3, 4], "term": [3, 4], "lie": 3, "group": 3, "g": [3, 4], "act": [3, 4], "transit": 3, "phase": 3, "space": [3, 4], "mathcal": 3, "present": 3, "infinitesim": 3, "gener": [3, 4], "action": [3, 4], "_eom": 3, "let": 3, "u": [3, 4], "consid": [3, 4], "subject": 3, "graviti": 3, "model": [3, 4], "rigid": [3, 4], "massless": 3, "link": 3, "serial": 3, "spheric": 3, "joint": [3, 4], "first": [3, 4], "fix": 3, "point": 3, "place": 3, "origin": 3, "ambient": 3, "mathbb": 3, "neglect": 3, "friction": 3, "interact": 3, "among": 3, "lee": 3, "leok": 3, "mcclamroch": 3, "2018": 3, "omit": 3, "detail": 3, "q_i": 3, "": [3, 4], "configur": [3, 4], "vector": [3, 4], "th": 3, "mass": [3, 4], "m_i": 3, "express": 3, "our": 3, "q_1": 3, "dot": [3, 4], "q_n": 3, "subset": 3, "3n": 3, "angular": 3, "veloc": [3, 4], "omega_1": 3, "omega_n": 3, "t_": [3, 4], "time": [3, 4], "defin": 3, "kinemat": [3, 4], "begin": [3, 4], "align": [3, 4], "q": [3, 4], "_i": [3, 4], "omega_i": 3, "quad": 3, "1": [3, 4], "end": [3, 4], "written": [3, 4], "where": [3, 4], "symmetr": [3, 4], "block": 3, "matrix": [3, 4], "_": [3, 4], "ii": 3, "big": 3, "sum_": [3, 4], "j": 3, "nm_j": 3, "l_i": 3, "2i_3": 3, "ij": 3, "k": [3, 4], "m_k": 3, "l_il_j": 3, "hat": [3, 4], "_j": 3, "ji": 3, "m_": 3, "text": 3, "max": 3, "i_3": 3, "henc": 3, "field": 3, "f": [3, 4], "mathfrak": 3, "x": 3, "now": 3, "find": 3, "rightarrow": 3, "psi_": 3, "vert_m": 3, "foral": 3, "psi": 3, "subsect": 3, "ref": 3, "286subsec": 3, "cartesian": 3, "sinc": [3, 4], "linear": [3, 4], "invert": 3, "map": 3, "a_": [3, 4], "a_q": 3, "omega": 3, "rewrit": 3, "od": 3, "left": [3, 4], "bmatrix": [3, 4], "r_1": 3, "vdot": 3, "r_n": 3, "right": [3, 4], "h_1": 3, "h_n": 3, "a_1": 3, "a_n": 3, "In": [3, 4], "r_i": 3, "a_i": 3, "h_i": 3, "thu": 3, "given": [3, 4], "simeq": 3, "6n": 3, "electromechan": 4, "coupl": 4, "beam": 4, "within": 4, "constrain": 4, "scheme": 4, "null": 4, "d": 4, "alembert": 4, "principl": 4, "extend": 4, "enforc": 4, "constraint": 4, "via": 4, "multipli": 4, "delta": 4, "int_": 4, "0": 4, "l": 4, "mathbf": 4, "cdot": 4, "boldsymbol": 4, "lambda": 4, "dt": 4, "rm": 4, "ext": 4, "lagrangian": 4, "repres": 4, "holonom": 4, "extern": 4, "forc": 4, "By": 4, "electr": 4, "effect": 4, "geometr": 4, "exact": 4, "potenti": 4, "phi_o": 4, "increment": 4, "alpha": 4, "beta": 4, "treat": 4, "degre": 4, "freedom": 4, "phi": 4, "varphi": 4, "_1": 4, "_2": 4, "_3": 4, "accord": 4, "assumpt": 4, "director": 4, "fulfil": 4, "orthogon": 4, "frac": 4, "differ": 4, "between": 4, "kinet": 4, "energi": 4, "v": 4, "do": 4, "contribut": 4, "int_c": 4, "rho": 4, "i_": 4, "label": 4, "densiti": 4, "per": 4, "arc": 4, "length": 4, "moment": 4, "inertia": 4, "cross": 4, "section": 4, "compon": 4, "consist": 4, "correspond": 4, "zero": 4, "hyperelast": 4, "materi": 4, "dea": 4, "an": 4, "strain": 4, "omega_b": 4, "over": 4, "center": 4, "line": 4, "all": 4, "non": 4, "conserv": 4, "viscoelast": 4, "base": 4, "kelvin": 4, "voigt": 4, "w": 4, "vi": 4, "b_0": 4, "p": 4, "dv": 4, "two": 4, "conjug": 4, "quantiti": 4, "being": 4, "piola": 4, "kirchhoff": 4, "stress": 4, "deform": 4, "gradient": 4, "case": 4, "formul": 4, "partial": 4, "sigma": 4, "dad": 4, "spatial": 4, "1d": 4, "finit": 4, "element": 4, "one": 4, "dimension": 4, "type": 4, "shape": 4, "appli": 4, "directli": 4, "togeth": 4, "centroid": 4, "Then": 4, "tempor": 4, "which": 4, "good": 4, "long": 4, "behavior": 4, "interv": 4, "t_n": 4, "l_d": 4, "approx": 4, "_n": 4, "midpoint": 4, "rule": 4, "mid": 4, "after": 4, "take": 4, "stationar": 4, "elimin": 4, "nodal": 4, "reparametr": 4, "_d": 4, "lead": 4, "unknown": 4, "evalu": 4, "int": 4, "skew": 4, "ident": 4, "multibodi": 4, "compos": 4, "flexibl": 4, "actuat": 4, "bodi": 4, "design": 4, "extra": 4, "well": 4, "solv": 4, "effici": 4, "reduc": 4, "further": 4, "minim": 4, "possibl": 4, "dimens": 4, "specifi": 4, "theta": 4, "character": 4, "displac": 4, "rotat": 4, "respect": 4, "next": 4, "updat": 4, "exp": 4, "mean": 4, "chang": 4, "nonlinear": 4, "newton": 4, "rapson": 4, "algorithm": 4, "tangent": 4, "iter": 4, "_t": 4, "residu": 4}, "objects": {"": [[2, 0, 0, "-", "src"], [2, 0, 1, "", "src"]], "src": [[2, 1, 1, "", "RKMK45"], [2, 1, 1, "", "RKMK5"], [2, 2, 1, "", "main_const_vs_var_stepsize_comparison"]]}, "objtypes": {"0": "mat:module", "1": "mat:script", "2": "mat:function"}, "objnames": {"0": ["mat", "module", "MATLAB module"], "1": ["mat", "script", "MATLAB script"], "2": ["mat", "function", "MATLAB function"]}, "titleterms": {"document": 0, "variable_stepsize_lie_group_integr": 0, "content": 0, "indic": 0, "tabl": 0, "instal": 1, "virtual": 1, "environ": 1, "matlab": 2, "code": 2, "src": 2, "The": 3, "n": 3, "fold": 3, "3d": 3, "pendulum": 3, "equat": [3, 4], "motion": 3, "variat": 4, "integr": 4, "continu": 4, "set": 4, "discret": 4, "euler": 4, "lagrang": 4}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx": 56}})