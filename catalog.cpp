#include "phz.hpp"

namespace phz {
const std::map<std::string, std::string>& examples() {
    static const std::map<std::string, std::string> d = {
        {"aq1", R"PHZMODEL({
  "nm": "Automatic quadrature of exponential reaction load",
  "dim": 1,
  "fld": [
    "u"
  ],
  "eq": {
    "u": "-u+exp(5*x)"
  },
  "ini": {
    "u": 0
  },
  "dom": {
    "type": "box",
    "n": [
      2
    ]
  },
  "disc": {
    "type": "fe",
    "q": 2,
    "aq": {
      "atol": 1e-08,
      "rtol": 1e-07
    }
  },
  "stat": {
    "tol": 1e-10,
    "eig": false
  },
  "run": {
    "method": "BE",
    "t": [
      0,
      0.1
    ],
    "n": 3,
    "max_step": 0.05,
    "rtol": 1e-10,
    "atol": 1e-12
  },
  "asm": [
    "Stationary state is the P1 L2 projection of exp(5*x)",
    "Quadrature control estimates sensitivity with mesh and time settings fixed"
  ]
})PHZMODEL"},
        {"aq3", R"PHZMODEL({
  "nm": "Automatic quadrature of exponential reaction load",
  "dim": 3,
  "fld": [
    "u"
  ],
  "eq": {
    "u": "-u+exp(5*x)"
  },
  "ini": {
    "u": 0
  },
  "dom": {
    "type": "box",
    "n": [
      2,
      2,
      2
    ]
  },
  "disc": {
    "type": "fe",
    "q": 3,
    "aq": {
      "atol": 1e-08,
      "rtol": 1e-07
    }
  },
  "stat": {
    "tol": 1e-10,
    "eig": false
  },
  "run": {
    "method": "BE",
    "t": [
      0,
      0.1
    ],
    "n": 3,
    "max_step": 0.05,
    "rtol": 1e-10,
    "atol": 1e-12
  },
  "asm": [
    "Stationary state is the P1 L2 projection of exp(5*x)",
    "Quadrature control estimates sensitivity with mesh and time settings fixed"
  ]
})PHZMODEL"},
        {"aqns2", R"PHZMODEL({
  "dim": 2,
  "nm": "Automatic quadrature of smooth nonlinear incompressible flow",
  "fld": {
    "u": [
      2
    ],
    "p": []
  },
  "par": {
    "nu": 0.5
  },
  "coef": {
    "f": [
      2
    ]
  },
  "src": {
    "f": "[12*pi**3*nu*sin(pi*x)**2*sin(pi*y)*cos(pi*y) - 4*pi**3*nu*sin(pi*y)*cos(pi*x)**2*cos(pi*y) + 4*pi**3*sin(pi*x)**3*sin(pi*y)**4*cos(pi*x) + 4*pi**3*sin(pi*x)**3*sin(pi*y)**2*cos(pi*x)*cos(pi*y)**2 - pi*sin(pi*x)*cos(pi*y),-12*pi**3*nu*sin(pi*x)*sin(pi*y)**2*cos(pi*x) + 4*pi**3*nu*sin(pi*x)*cos(pi*x)*cos(pi*y)**2 + 4*pi**3*sin(pi*x)**4*sin(pi*y)**3*cos(pi*y) + 4*pi**3*sin(pi*x)**2*sin(pi*y)**3*cos(pi*x)**2*cos(pi*y) - pi*sin(pi*y)*cos(pi*x)]"
  },
  "pde": {
    "u": "dt(u)+adv(u,u)=-grad(p)+nu*lap(u)+f",
    "p": "div(u)=0"
  },
  "ini": {
    "u": 0,
    "p": 0
  },
  "bc": {
    "u": {
      "type": "dirichlet",
      "val": "[2*pi*sin(pi*x)**2*sin(pi*y)*cos(pi*y),-2*pi*sin(pi*x)*sin(pi*y)**2*cos(pi*x)]"
    }
  },
  "dom": {
    "type": "box",
    "n": [
      3,
      3
    ]
  },
  "disc": {
    "type": "th",
    "q": 4,
    "aq": {
      "atol": 1e-05,
      "rtol": 0.0001,
      "step": 1
    }
  },
  "run": {
    "method": "BE",
    "t": [
      0,
      0.1
    ],
    "n": 3,
    "max_step": 0.05
  },
  "stat": {
    "tol": 1e-10
  },
  "chk": {
    "exact": {
      "u": "[2*pi*sin(pi*x)**2*sin(pi*y)*cos(pi*y),-2*pi*sin(pi*x)*sin(pi*y)**2*cos(pi*x)]",
      "p": "cos(pi*x)*cos(pi*y)"
    },
    "tol": 1.0
  },
  "asm": [
    "Manufactured incompressible solution; source generated from the strong PDE",
    "Continuous P2 velocity and P1 pressure; pressure has zero volume mean",
    "Full Dirichlet velocity; nonlinear terms use element quadrature",
    "Use Equilibrium for stationary examples; run for transient examples"
  ]
}
)PHZMODEL"},
        {"aqns3", R"PHZMODEL({
  "dim": 3,
  "nm": "Automatic quadrature of smooth nonlinear incompressible flow",
  "fld": {
    "u": [
      3
    ],
    "p": []
  },
  "par": {
    "nu": 0.5
  },
  "coef": {
    "f": [
      3
    ]
  },
  "src": {
    "f": "[16*pi**3*nu*sin(pi*x)**2*sin(pi*y)*sin(pi*z)**2*cos(pi*y) - 4*pi**3*nu*sin(pi*x)**2*sin(pi*y)*cos(pi*y)*cos(pi*z)**2 - 4*pi**3*nu*sin(pi*y)*sin(pi*z)**2*cos(pi*x)**2*cos(pi*y) + 4*pi**3*sin(pi*x)**3*sin(pi*y)**4*sin(pi*z)**4*cos(pi*x) + 4*pi**3*sin(pi*x)**3*sin(pi*y)**2*sin(pi*z)**4*cos(pi*x)*cos(pi*y)**2 - pi*sin(pi*x)*cos(pi*y)*cos(pi*z),-16*pi**3*nu*sin(pi*x)*sin(pi*y)**2*sin(pi*z)**2*cos(pi*x) + 4*pi**3*nu*sin(pi*x)*sin(pi*y)**2*cos(pi*x)*cos(pi*z)**2 + 4*pi**3*nu*sin(pi*x)*sin(pi*z)**2*cos(pi*x)*cos(pi*y)**2 + 4*pi**3*sin(pi*x)**4*sin(pi*y)**3*sin(pi*z)**4*cos(pi*y) + 4*pi**3*sin(pi*x)**2*sin(pi*y)**3*sin(pi*z)**4*cos(pi*x)**2*cos(pi*y) - pi*sin(pi*y)*cos(pi*x)*cos(pi*z),-pi*sin(pi*z)*cos(pi*x)*cos(pi*y)]"
  },
  "pde": {
    "u": "dt(u)+adv(u,u)=-grad(p)+nu*lap(u)+f",
    "p": "div(u)=0"
  },
  "ini": {
    "u": 0,
    "p": 0
  },
  "bc": {
    "u": {
      "type": "dirichlet",
      "val": "[2*pi*sin(pi*x)**2*sin(pi*y)*sin(pi*z)**2*cos(pi*y),-2*pi*sin(pi*x)*sin(pi*y)**2*sin(pi*z)**2*cos(pi*x),0]"
    }
  },
  "dom": {
    "type": "box",
    "n": [
      2,
      2,
      2
    ]
  },
  "disc": {
    "type": "th",
    "q": 4,
    "aq": {
      "atol": 1e-05,
      "rtol": 0.0001,
      "step": 1,
      "qmax": 12
    }
  },
  "run": {
    "method": "BE",
    "t": [
      0,
      0.1
    ],
    "n": 3,
    "max_step": 0.05
  },
  "stat": {
    "tol": 1e-10
  },
  "chk": {
    "exact": {
      "u": "[2*pi*sin(pi*x)**2*sin(pi*y)*sin(pi*z)**2*cos(pi*y),-2*pi*sin(pi*x)*sin(pi*y)**2*sin(pi*z)**2*cos(pi*x),0]",
      "p": "cos(pi*x)*cos(pi*y)*cos(pi*z)"
    },
    "tol": 4.0
  },
  "asm": [
    "Manufactured smooth solution for mesh sensitivity",
    "Use a stationary mesh study",
    "Weak divergence moments and pointwise divergence are reported separately"
  ]
}
)PHZMODEL"},
        {"aqsp", R"PHZMODEL({
  "nm": "Automatic spectral quadrature of exponential load",
  "dim": 1,
  "fld": [
    "u"
  ],
  "eq": {
    "u": "-u+exp(5*x)"
  },
  "ini": {
    "u": 0
  },
  "dom": {
    "type": "box",
    "n": [
      2
    ]
  },
  "disc": {
    "type": "sp",
    "n": [
      2
    ],
    "q": 4,
    "aq": {
      "atol": 1e-09,
      "rtol": 1e-08,
      "qmax": 20,
      "step": 2
    }
  },
  "stat": {
    "tol": 1e-10,
    "eig": false
  },
  "run": {
    "method": "BE",
    "t": [
      0,
      0.1
    ],
    "n": 3,
    "max_step": 0.05,
    "rtol": 1e-10,
    "atol": 1e-12
  },
  "asm": [
    "Stationary state is the P1 L2 projection of exp(5*x)",
    "Quadrature control estimates sensitivity with mesh and time settings fixed"
  ]
})PHZMODEL"},
        {"be", R"PHZMODEL({
  "ver": 1,
  "nm": "Heat equation: explicit Newton residual history",
  "dim": 1,
  "fld": ["u"],
  "par": {"D": 0.2},
  "eq": {"u": "D*lap(u)"},
  "ini": {"u": "1+0.2*cos(pi*x)"},
  "dom": {"type": "box"},
  "disc": {"type": "sp", "n": [3]},
  "run": {"t": [0, 0.2], "n": 11, "method": "BE", "max_step": 0.01, "rtol": 1e-9, "atol": 1e-11},
  "chk": {"mass": ["u"], "exact": {"u": "1+0.2*exp(-D*pi**2*t)*cos(pi*x)"}, "tol": 0.001},
  "asm": ["Backward Euler is first order in time", "Every Newton iteration and every accepted step is saved", "rtol and atol control Newton residual, not the temporal truncation error"]
}
)PHZMODEL"},
        {"ch", R"PHZMODEL({
  "ver": 1,
  "nm": "Cahn-Hilliard: variational phase separation",
  "dim": 2,
  "fld": ["u"],
  "par": {"kap": 0.01, "M": 0.05},
  "en": {"f": "(u**2-1)**2/4", "kap": {"u": "kap"}, "mob": {"u": "M"}, "kin": {"u": "ch"}},
  "ini": {"u": "0.1+0.03*cos(pi*x)*cos(pi*y)"},
  "dom": {"type": "box"},
  "disc": {"type": "sp", "n": [3, 3], "q": 20},
  "run": {"t": [0, 2], "n": 31, "method": "BDF", "rtol": 1e-8, "atol": 1e-10},
  "chk": {"mass": ["u"], "energy": true, "tol": 1e-7},
  "asm": ["Dimensionless quartic free energy", "Constant mobility", "Natural zero-flux conditions for field and chemical potential", "A low-mode demonstration, not a converged interface study"]
}
)PHZMODEL"},
        {"dsl2", R"PHZMODEL({
  "ver": 1,
  "nm": "Nonlinear diffusion from an ASCII PDE with mixed boundaries",
  "dim": 2,
  "fld": [
    "u"
  ],
  "par": {
    "D": 0.2
  },
  "ini": {
    "u": "1+x+2*y"
  },
  "bc": {
    "u": [
      {
        "type": "dirichlet",
        "on": "x0",
        "val": "1+2*y+t**2"
      },
      {
        "type": "neumann",
        "on": "x1",
        "val": "-D*(1+(1+x+2*y+t**2)**2)"
      },
      {
        "type": "robin",
        "on": "y0",
        "coef": 2,
        "val": "(1+x+2*y+t**2)-D*(1+(1+x+2*y+t**2)**2)"
      },
      {
        "type": "neumann",
        "on": "y1",
        "val": "-2*D*(1+(1+x+2*y+t**2)**2)"
      }
    ]
  },
  "dom": {
    "type": "box",
    "n": [
      6,
      6
    ]
  },
  "disc": {
    "type": "fe",
    "q": 4
  },
  "run": {
    "t": [
      0,
      0.2
    ],
    "n": 11,
    "method": "BDF",
    "rtol": 1e-09,
    "atol": 1e-11
  },
  "chk": {
    "exact": {
      "u": "1+x+2*y+t**2"
    },
    "tol": 1e-07
  },
  "asm": [
    "Manufactured nonlinear diffusion with affine spatial solution",
    "D*(1+u**2) is evaluated at quadrature points before assembly",
    "Flux values include the complete nonlinear diffusion coefficient"
  ],
  "pde": {
    "u": "dt(u)-div(D*(1+u**2)*grad(u))=2*t-10*D*(1+x+2*y+t**2)"
  }
}
)PHZMODEL"},
        {"eq2", R"PHZMODEL({
  "ver": 1,
  "nm": "Nonlinear equilibrium with mixed boundaries in 2D",
  "dim": 2,
  "fld": ["u"],
  "par": {"D": 0.2},
  "eq": {"u": "D*lap(u)-u**3+(1+x+2*y)**3"},
  "ini": {"u": "1"},
  "bc": {"u": [
    {"type": "dirichlet", "on": "x0", "val": "1+2*y"},
    {"type": "neumann", "on": "x1", "val": "-D"},
    {"type": "robin", "on": "y0", "coef": 2, "val": "1+x-D"},
    {"type": "neumann", "on": "y1", "val": "-2*D"}
  ]},
  "dom": {"type": "box", "n": [6, 6]},
  "disc": {"type": "fe", "q": 4},
  "stat": {"tol": 1e-9, "max_it": 40, "eig": true},
  "chk": {"exact": {"u": "1+x+2*y"}, "nonneg": ["u"], "tol": 1e-7},
  "asm": ["Manufactured nonlinear equilibrium", "Initial field is a Newton starting guess", "Local stability is assessed on free FEM coefficients"]
}
)PHZMODEL"},
        {"eq3", R"PHZMODEL({
  "ver": 1,
  "nm": "Nonlinear equilibrium with mixed boundaries in 3D",
  "dim": 3,
  "fld": ["u"],
  "par": {"D": 0.2},
  "eq": {"u": "D*lap(u)-u**3+(1+x+2*y+3*z)**3"},
  "ini": {"u": "1"},
  "bc": {"u": [
    {"type": "dirichlet", "on": "x0", "val": "1+2*y+3*z"},
    {"type": "neumann", "on": "x1", "val": "-D"},
    {"type": "robin", "on": "y0", "coef": 2, "val": "1+x+3*z-D"},
    {"type": "neumann", "on": "y1", "val": "-2*D"},
    {"type": "neumann", "on": "z0", "val": "3*D"},
    {"type": "neumann", "on": "z1", "val": "-3*D"}
  ]},
  "dom": {"type": "box", "n": [3, 3, 3]},
  "disc": {"type": "fe", "q": 4},
  "stat": {"tol": 1e-9, "max_it": 40, "eig": true},
  "chk": {"exact": {"u": "1+x+2*y+3*z"}, "nonneg": ["u"], "tol": 1e-7},
  "asm": ["Manufactured nonlinear equilibrium on tetrahedra", "Boundary flux is positive outward", "Stability refers to the semidiscrete model at fixed boundary values"]
}
)PHZMODEL"},
        {"eqp", R"PHZMODEL({
  "ver": 1,
  "nm": "Closed electrolyte equilibrium at fixed ion masses",
  "dim": 2,
  "fld": ["cp", "cm"],
  "alg": [{"nm": "psi", "a": "eps", "src": "cp-cm"}],
  "par": {"D": 0.1, "eps": 0.5},
  "eq": {"cp": "D*lap(cp)+D*divg(cp,psi)", "cm": "D*lap(cm)-D*divg(cm,psi)"},
  "ini": {"cp": "1+0.05*cos(pi*x)", "cm": "1-0.05*cos(pi*x)"},
  "dom": {"type": "box", "n": [3, 3]},
  "disc": {"type": "fe", "q": 4},
  "stat": {"mass": {"cp": "initial", "cm": "initial"}, "tol": 1e-9, "eig": true},
  "chk": {"exact": {"cp": "1", "cm": "1"}, "nonneg": ["cp", "cm"], "tol": 1e-7},
  "asm": ["Dimensionless monovalent electrolyte without fixed charge", "Both ion masses are fixed to their starting values", "All Newton perturbations preserve both masses and charge neutrality", "Potential has zero spatial mean"]
}
)PHZMODEL"},
        {"eqs2", R"PHZMODEL({
  "ver": 1,
  "nm": "Sparse nonlinear equilibrium mesh study in 2D",
  "dim": 2,
  "fld": ["u"],
  "eq": {"u": "lap(u)-u**3+2*pi**2*sin(pi*x)*sin(pi*y)+(sin(pi*x)*sin(pi*y))**3"},
  "ini": {"u": 0},
  "bc": {"u": {"type": "dirichlet", "val": 0}},
  "dom": {"type": "box", "n": [5, 5]},
  "disc": {"type": "fe", "q": 4},
  "stat": {"solver": "sparse", "tol": 1e-9, "eig": false},
  "chk": {"exact": {"u": "sin(pi*x)*sin(pi*y)"}, "tol": 0.1},
  "asm": ["Manufactured smooth equilibrium for independent mesh convergence", "Use a three-level stationary mesh study; quadrature remains fixed", "The coarse exact-solution tolerance is not a solver tolerance"]
}
)PHZMODEL"},
        {"eqs3", R"PHZMODEL({
  "ver": 1,
  "nm": "Sparse nonlinear diffusion equilibrium in 3D",
  "dim": 3,
  "fld": ["u"],
  "eq": {"u": "divg(1+u**2,u)-u**3+(1+x+y+z)**3-6*(1+x+y+z)"},
  "ini": {"u": "1"},
  "bc": {"u": {"type": "dirichlet", "val": "1+x+y+z"}},
  "dom": {"type": "box", "n": [9, 9, 9]},
  "disc": {"type": "fe", "q": 4},
  "stat": {"solver": "sparse", "tol": 1e-8, "max_it": 40, "eig": true, "timeout": 120},
  "chk": {"exact": {"u": "1+x+y+z"}, "nonneg": ["u"], "tol": 1e-7},
  "asm": ["Manufactured nonlinear diffusion and reaction on tetrahedra", "The exact affine field belongs to the P1 space", "A large equilibrium does not by itself establish stability; spectrum above 256 free coefficients is not computed"]
}
)PHZMODEL"},
        {"fh", R"PHZMODEL({
  "ver": 1,
  "nm": "Flory-Huggins mixing: conserved dynamics",
  "dim": 2,
  "fld": ["phi"],
  "par": {"Np": 10, "Ns": 1, "chi": 1.0, "kap": 0.005, "M": 0.01},
  "en": {"f": "phi*log(phi)/Np+(1-phi)*log(1-phi)/Ns+chi*phi*(1-phi)", "kap": {"phi": "kap"}, "mob": {"phi": "M"}, "kin": {"phi": "ch"}},
  "ini": {"phi": "0.4+0.01*cos(pi*x)*cos(pi*y)"},
  "dom": {"type": "box", "n": [6, 6]},
  "disc": {"type": "fe", "q": 4},
  "run": {"t": [0, 0.2], "n": 11, "method": "BDF", "rtol": 1e-7, "atol": 1e-9},
  "chk": {"mass": ["phi"], "energy": true, "nonneg": ["phi"], "tol": 1e-6},
  "asm": ["Mixing and gradient energy only", "0 < phi < 1 is the logarithm domain", "No elastic deformation or ion coupling in this example", "Mixed P1 discretization of chemical potential"]
}
)PHZMODEL"},
        {"heat3", R"PHZMODEL({
  "ver": 1,
  "nm": "3D heat equation: exact cosine mode",
  "dim": 3,
  "fld": ["u"],
  "par": {"D": 0.1},
  "eq": {"u": "D*lap(u)"},
  "ini": {"u": "1+0.2*cos(pi*x)*cos(pi*y)*cos(pi*z)"},
  "dom": {"type": "box", "n": [3, 3, 3]},
  "disc": {"type": "fe", "q": 4},
  "run": {"t": [0, 0.2], "n": 11, "method": "BDF", "rtol": 1e-8, "atol": 1e-10},
  "chk": {"mass": ["u"], "nonneg": ["u"], "exact": {"u": "1+0.2*exp(-3*pi**2*D*t)*cos(pi*x)*cos(pi*y)*cos(pi*z)"}, "tol": 0.04},
  "asm": ["Unit cube", "Homogeneous Neumann boundary", "Continuous linear tetrahedra", "Exact solution is used to measure spatial and temporal error"]
}
)PHZMODEL"},
        {"int3", R"PHZMODEL({
  "ver": 1,
  "nm": "3D nonlocal diffusion with conserved mass",
  "dim": 3,
  "fld": [
    "u"
  ],
  "par": {
    "D": 0.1,
    "b": 0.3
  },
  "ini": {
    "u": "1+0.2*cos(pi*x)*cos(pi*y)*cos(pi*z)"
  },
  "dom": {
    "type": "box",
    "n": [
      3,
      3,
      3
    ]
  },
  "disc": {
    "type": "fe",
    "q": 4
  },
  "run": {
    "t": [
      0,
      0.2
    ],
    "n": 11,
    "method": "BDF",
    "rtol": 1e-08,
    "atol": 1e-10
  },
  "chk": {
    "mass": [
      "u"
    ],
    "nonneg": [
      "u"
    ],
    "exact": {
      "u": "1+0.2*exp(-(3*pi**2*D+b)*t)*cos(pi*x)*cos(pi*y)*cos(pi*z)"
    },
    "tol": 0.04
  },
  "asm": [
    "Unit cube with zero outward flux",
    "intg is evaluated by volume quadrature on the actual mesh",
    "The nonlocal exchange term conserves total mass",
    "Mesh convergence and quadrature sensitivity are independent checks"
  ],
  "pde": {
    "u": "dt(u)=D*lap(u)+b*(intg(u)/intg(1)-u)"
  }
}
)PHZMODEL"},
        {"lshape", R"PHZMODEL({
  "ver": 1,
  "nm": "Reaction-diffusion on an L-shaped domain",
  "dim": 2,
  "fld": ["u", "v"],
  "par": {"Du": 0.02, "Dv": 0.08, "a": 1.0, "b": 0.5},
  "eq": {"u": "Du*lap(u)+a*u-b*u*v", "v": "Dv*lap(v)+b*u*v-v"},
  "ini": {"u": "2+0.3*exp(-40*((x-0.25)**2+(y-0.7)**2))", "v": "1+0.1*cos(pi*x)"},
  "dom": {"type": "lshape", "n": [8, 8]},
  "disc": {"type": "fe", "q": 4},
  "run": {"t": [0, 2], "n": 31, "method": "BDF", "rtol": 1e-7, "atol": 1e-9},
  "chk": {"nonneg": ["u", "v"], "tol": 1e-6},
  "asm": ["Fixed nonconvex polygon", "Continuous P1 triangles", "No boundary flux", "Reaction terms integrated at element quadrature points"]
}
)PHZMODEL"},
        {"mix2", R"PHZMODEL({
  "ver": 1,
  "nm": "Mixed boundaries: time-dependent Dirichlet, flux and Robin",
  "dim": 2,
  "fld": ["u"],
  "par": {"D": 0.2},
  "eq": {"u": "D*lap(u)+2*t"},
  "ini": {"u": "1+x+2*y"},
  "bc": {"u": [
    {"type": "dirichlet", "on": "x0", "val": "1+2*y+t**2"},
    {"type": "neumann", "on": "x1", "val": "-D"},
    {"type": "robin", "on": "y0", "coef": 2, "val": "1+x+t**2-D"},
    {"type": "neumann", "on": "y1", "val": "-2*D"}
  ]},
  "dom": {"type": "box", "n": [6, 6]},
  "disc": {"type": "fe", "q": 4},
  "run": {"t": [0, 0.2], "n": 11, "method": "BDF", "rtol": 1e-9, "atol": 1e-11},
  "chk": {"exact": {"u": "1+x+2*y+t**2"}, "tol": 1e-7},
  "asm": ["Manufactured solution is affine in space and quadratic in time", "Positive prescribed flux is outward and already includes D", "Dirichlet values are lifted exactly; free equations include -Mfc*dg/dt"]
}
)PHZMODEL"},
        {"mix3", R"PHZMODEL({
  "ver": 1,
  "nm": "3D mixed boundaries on a tetrahedral mesh",
  "dim": 3,
  "fld": ["u"],
  "par": {"D": 0.2},
  "eq": {"u": "D*lap(u)+2*t"},
  "ini": {"u": "1+x+2*y+3*z"},
  "bc": {"u": [
    {"type": "dirichlet", "on": "x0", "val": "1+2*y+3*z+t**2"},
    {"type": "neumann", "on": "x1", "val": "-D"},
    {"type": "robin", "on": "y0", "coef": 2, "val": "1+x+3*z+t**2-D"},
    {"type": "neumann", "on": "y1", "val": "-2*D"},
    {"type": "neumann", "on": "z0", "val": "3*D"},
    {"type": "neumann", "on": "z1", "val": "-3*D"}
  ]},
  "dom": {"type": "box", "n": [3, 3, 3]},
  "disc": {"type": "fe", "q": 4},
  "run": {"t": [0, 0.2], "n": 11, "method": "BDF", "rtol": 1e-9, "atol": 1e-11},
  "chk": {"exact": {"u": "1+x+2*y+3*z+t**2"}, "tol": 1e-7},
  "asm": ["Manufactured affine spatial solution on P1 tetrahedra", "Natural loads are integrated on triangular boundary facets", "Unspecified boundary facets carry zero total outward flux"]
}
)PHZMODEL"},
        {"nseq2", R"PHZMODEL({
  "dim": 2,
  "nm": "Stationary nonlinear incompressible flow in 2D",
  "fld": {
    "u": [
      2
    ],
    "p": []
  },
  "par": {
    "nu": 0.5
  },
  "coef": {
    "f": [
      2
    ]
  },
  "src": {
    "f": "[-2*nu + 2*x**3 + 1,2*x**2*y + 1]"
  },
  "pde": {
    "u": "dt(u)+adv(u,u)=-grad(p)+nu*lap(u)+f",
    "p": "div(u)=0"
  },
  "ini": {
    "u": 0,
    "p": 0
  },
  "bc": {
    "u": {
      "type": "dirichlet",
      "val": "[x**2,-2*x*y]"
    }
  },
  "dom": {
    "type": "box",
    "n": [
      3,
      3
    ]
  },
  "disc": {
    "type": "th",
    "q": 4
  },
  "run": {
    "method": "BE",
    "t": [
      0,
      0.1
    ],
    "n": 3,
    "max_step": 0.05
  },
  "stat": {
    "tol": 1e-10
  },
  "chk": {
    "exact": {
      "u": "[x**2,-2*x*y]",
      "p": "x + y - 1"
    },
    "tol": 1e-07
  },
  "asm": [
    "Manufactured incompressible solution; source generated from the strong PDE",
    "Continuous P2 velocity and P1 pressure; pressure has zero volume mean",
    "Full Dirichlet velocity; nonlinear terms use element quadrature",
    "Use Equilibrium for stationary examples; run for transient examples"
  ]
}
)PHZMODEL"},
        {"nseq3", R"PHZMODEL({
  "dim": 3,
  "nm": "Stationary nonlinear incompressible flow in 3D",
  "fld": {
    "u": [
      3
    ],
    "p": []
  },
  "par": {
    "nu": 0.5
  },
  "coef": {
    "f": [
      3
    ]
  },
  "src": {
    "f": "[-2*nu + 2*x**3 + 1,1,1]"
  },
  "pde": {
    "u": "dt(u)+adv(u,u)=-grad(p)+nu*lap(u)+f",
    "p": "div(u)=0"
  },
  "ini": {
    "u": 0,
    "p": 0
  },
  "bc": {
    "u": {
      "type": "dirichlet",
      "val": "[x**2,-x*y,-x*z]"
    }
  },
  "dom": {
    "type": "box",
    "n": [
      2,
      2,
      2
    ]
  },
  "disc": {
    "type": "th",
    "q": 4
  },
  "run": {
    "method": "BE",
    "t": [
      0,
      0.1
    ],
    "n": 3,
    "max_step": 0.05
  },
  "stat": {
    "tol": 1e-10
  },
  "chk": {
    "exact": {
      "u": "[x**2,-x*y,-x*z]",
      "p": "x + y + z - 3/2"
    },
    "tol": 1e-07
  },
  "asm": [
    "Manufactured incompressible solution; source generated from the strong PDE",
    "Continuous P2 velocity and P1 pressure; pressure has zero volume mean",
    "Full Dirichlet velocity; nonlinear terms use element quadrature",
    "Use Equilibrium for stationary examples; run for transient examples"
  ]
}
)PHZMODEL"},
        {"nsref2", R"PHZMODEL({
  "dim": 2,
  "nm": "Stationary nonlinear incompressible flow in 2D",
  "fld": {
    "u": [
      2
    ],
    "p": []
  },
  "par": {
    "nu": 0.5
  },
  "coef": {
    "f": [
      2
    ]
  },
  "src": {
    "f": "[12*pi**3*nu*sin(pi*x)**2*sin(pi*y)*cos(pi*y) - 4*pi**3*nu*sin(pi*y)*cos(pi*x)**2*cos(pi*y) + 4*pi**3*sin(pi*x)**3*sin(pi*y)**4*cos(pi*x) + 4*pi**3*sin(pi*x)**3*sin(pi*y)**2*cos(pi*x)*cos(pi*y)**2 - pi*sin(pi*x)*cos(pi*y),-12*pi**3*nu*sin(pi*x)*sin(pi*y)**2*cos(pi*x) + 4*pi**3*nu*sin(pi*x)*cos(pi*x)*cos(pi*y)**2 + 4*pi**3*sin(pi*x)**4*sin(pi*y)**3*cos(pi*y) + 4*pi**3*sin(pi*x)**2*sin(pi*y)**3*cos(pi*x)**2*cos(pi*y) - pi*sin(pi*y)*cos(pi*x)]"
  },
  "pde": {
    "u": "dt(u)+adv(u,u)=-grad(p)+nu*lap(u)+f",
    "p": "div(u)=0"
  },
  "ini": {
    "u": 0,
    "p": 0
  },
  "bc": {
    "u": {
      "type": "dirichlet",
      "val": "[2*pi*sin(pi*x)**2*sin(pi*y)*cos(pi*y),-2*pi*sin(pi*x)*sin(pi*y)**2*cos(pi*x)]"
    }
  },
  "dom": {
    "type": "box",
    "n": [
      3,
      3
    ]
  },
  "disc": {
    "type": "th",
    "q": 4
  },
  "run": {
    "method": "BE",
    "t": [
      0,
      0.1
    ],
    "n": 3,
    "max_step": 0.05
  },
  "stat": {
    "tol": 1e-10
  },
  "chk": {
    "exact": {
      "u": "[2*pi*sin(pi*x)**2*sin(pi*y)*cos(pi*y),-2*pi*sin(pi*x)*sin(pi*y)**2*cos(pi*x)]",
      "p": "cos(pi*x)*cos(pi*y)"
    },
    "tol": 1.0
  },
  "asm": [
    "Manufactured incompressible solution; source generated from the strong PDE",
    "Continuous P2 velocity and P1 pressure; pressure has zero volume mean",
    "Full Dirichlet velocity; nonlinear terms use element quadrature",
    "Use Equilibrium for stationary examples; run for transient examples"
  ]
}
)PHZMODEL"},
        {"nsref3", R"PHZMODEL({
  "dim": 3,
  "nm": "3D smooth incompressible flow: mesh study",
  "fld": {
    "u": [
      3
    ],
    "p": []
  },
  "par": {
    "nu": 0.5
  },
  "coef": {
    "f": [
      3
    ]
  },
  "src": {
    "f": "[16*pi**3*nu*sin(pi*x)**2*sin(pi*y)*sin(pi*z)**2*cos(pi*y) - 4*pi**3*nu*sin(pi*x)**2*sin(pi*y)*cos(pi*y)*cos(pi*z)**2 - 4*pi**3*nu*sin(pi*y)*sin(pi*z)**2*cos(pi*x)**2*cos(pi*y) + 4*pi**3*sin(pi*x)**3*sin(pi*y)**4*sin(pi*z)**4*cos(pi*x) + 4*pi**3*sin(pi*x)**3*sin(pi*y)**2*sin(pi*z)**4*cos(pi*x)*cos(pi*y)**2 - pi*sin(pi*x)*cos(pi*y)*cos(pi*z),-16*pi**3*nu*sin(pi*x)*sin(pi*y)**2*sin(pi*z)**2*cos(pi*x) + 4*pi**3*nu*sin(pi*x)*sin(pi*y)**2*cos(pi*x)*cos(pi*z)**2 + 4*pi**3*nu*sin(pi*x)*sin(pi*z)**2*cos(pi*x)*cos(pi*y)**2 + 4*pi**3*sin(pi*x)**4*sin(pi*y)**3*sin(pi*z)**4*cos(pi*y) + 4*pi**3*sin(pi*x)**2*sin(pi*y)**3*sin(pi*z)**4*cos(pi*x)**2*cos(pi*y) - pi*sin(pi*y)*cos(pi*x)*cos(pi*z),-pi*sin(pi*z)*cos(pi*x)*cos(pi*y)]"
  },
  "pde": {
    "u": "dt(u)+adv(u,u)=-grad(p)+nu*lap(u)+f",
    "p": "div(u)=0"
  },
  "ini": {
    "u": 0,
    "p": 0
  },
  "bc": {
    "u": {
      "type": "dirichlet",
      "val": "[2*pi*sin(pi*x)**2*sin(pi*y)*sin(pi*z)**2*cos(pi*y),-2*pi*sin(pi*x)*sin(pi*y)**2*sin(pi*z)**2*cos(pi*x),0]"
    }
  },
  "dom": {
    "type": "box",
    "n": [
      2,
      2,
      2
    ]
  },
  "disc": {
    "type": "th",
    "q": 4
  },
  "run": {
    "method": "BE",
    "t": [
      0,
      0.1
    ],
    "n": 3,
    "max_step": 0.05
  },
  "stat": {
    "tol": 1e-10
  },
  "chk": {
    "exact": {
      "u": "[2*pi*sin(pi*x)**2*sin(pi*y)*sin(pi*z)**2*cos(pi*y),-2*pi*sin(pi*x)*sin(pi*y)**2*sin(pi*z)**2*cos(pi*x),0]",
      "p": "cos(pi*x)*cos(pi*y)*cos(pi*z)"
    },
    "tol": 4.0
  },
  "asm": [
    "Manufactured smooth solution for mesh sensitivity",
    "Use a stationary mesh study",
    "Weak divergence moments and pointwise divergence are reported separately"
  ]
}
)PHZMODEL"},
        {"nst2", R"PHZMODEL({
  "dim": 2,
  "nm": "Transient nonlinear incompressible flow in 2D",
  "fld": {
    "u": [
      2
    ],
    "p": []
  },
  "par": {
    "nu": 0.5
  },
  "coef": {
    "f": [
      2
    ]
  },
  "src": {
    "f": "[-2*nu*t - 2*nu + 2*t**2*x**3 + 4*t*x**3 + t + 2*x**3 + x**2 + 1,2*t**2*x**2*y + 4*t*x**2*y + t + 2*x**2*y - 2*x*y + 1]"
  },
  "pde": {
    "u": "dt(u)+adv(u,u)=-grad(p)+nu*lap(u)+f",
    "p": "div(u)=0"
  },
  "ini": {
    "u": "[x**2,-2*x*y]",
    "p": "x + y - 1"
  },
  "bc": {
    "u": {
      "type": "dirichlet",
      "val": "[x**2*(t + 1),-2*x*y*(t + 1)]"
    }
  },
  "dom": {
    "type": "box",
    "n": [
      3,
      3
    ]
  },
  "disc": {
    "type": "th",
    "q": 4
  },
  "run": {
    "method": "BE",
    "t": [
      0,
      0.1
    ],
    "n": 3,
    "max_step": 0.05
  },
  "stat": {
    "tol": 1e-10
  },
  "chk": {
    "exact": {
      "u": "[x**2*(t + 1),-2*x*y*(t + 1)]",
      "p": "(t + 1)*(x + y - 1)"
    },
    "tol": 1e-07
  },
  "asm": [
    "Manufactured incompressible solution; source generated from the strong PDE",
    "Continuous P2 velocity and P1 pressure; pressure has zero volume mean",
    "Full Dirichlet velocity; nonlinear terms use element quadrature",
    "Use Equilibrium for stationary examples; run for transient examples"
  ]
}
)PHZMODEL"},
        {"nst3", R"PHZMODEL({
  "dim": 3,
  "nm": "Transient nonlinear incompressible flow in 3D",
  "fld": {
    "u": [
      3
    ],
    "p": []
  },
  "par": {
    "nu": 0.5
  },
  "coef": {
    "f": [
      3
    ]
  },
  "src": {
    "f": "[-2*nu*t - 2*nu + 2*t**2*x**3 + 4*t*x**3 + t + 2*x**3 + x**2 + 1,t - x*y + 1,t - x*z + 1]"
  },
  "pde": {
    "u": "dt(u)+adv(u,u)=-grad(p)+nu*lap(u)+f",
    "p": "div(u)=0"
  },
  "ini": {
    "u": "[x**2,-x*y,-x*z]",
    "p": "x + y + z - 3/2"
  },
  "bc": {
    "u": {
      "type": "dirichlet",
      "val": "[x**2*(t + 1),-x*y*(t + 1),-x*z*(t + 1)]"
    }
  },
  "dom": {
    "type": "box",
    "n": [
      2,
      2,
      2
    ]
  },
  "disc": {
    "type": "th",
    "q": 4
  },
  "run": {
    "method": "BE",
    "t": [
      0,
      0.1
    ],
    "n": 3,
    "max_step": 0.05
  },
  "stat": {
    "tol": 1e-10
  },
  "chk": {
    "exact": {
      "u": "[x**2*(t + 1),-x*y*(t + 1),-x*z*(t + 1)]",
      "p": "(t + 1)*(x + y + z - 3/2)"
    },
    "tol": 1e-07
  },
  "asm": [
    "Manufactured incompressible solution; source generated from the strong PDE",
    "Continuous P2 velocity and P1 pressure; pressure has zero volume mean",
    "Full Dirichlet velocity; nonlinear terms use element quadrature",
    "Use Equilibrium for stationary examples; run for transient examples"
  ]
}
)PHZMODEL"},
        {"nstime", R"PHZMODEL({
  "dim": 2,
  "nm": "Transient nonlinear incompressible flow in 2D",
  "fld": {
    "u": [
      2
    ],
    "p": []
  },
  "par": {
    "nu": 0.5
  },
  "coef": {
    "f": [
      2
    ]
  },
  "src": {
    "f": "[-2*nu*t**2 - 4*nu*t - 2*nu + 2*t**4*x**3 + 8*t**3*x**3 + 12*t**2*x**3 + t**2 + 8*t*x**3 + 2*t*x**2 + 2*t + 2*x**3 + 2*x**2 + 1,2*t**4*x**2*y + 8*t**3*x**2*y + 12*t**2*x**2*y + t**2 + 8*t*x**2*y - 4*t*x*y + 2*t + 2*x**2*y - 4*x*y + 1]"
  },
  "pde": {
    "u": "dt(u)+adv(u,u)=-grad(p)+nu*lap(u)+f",
    "p": "div(u)=0"
  },
  "ini": {
    "u": "[x**2,-2*x*y]",
    "p": "x + y - 1"
  },
  "bc": {
    "u": {
      "type": "dirichlet",
      "val": "[x**2*(t + 1)**2,-2*x*y*(t + 1)**2]"
    }
  },
  "dom": {
    "type": "box",
    "n": [
      3,
      3
    ]
  },
  "disc": {
    "type": "th",
    "q": 4
  },
  "run": {
    "method": "BE",
    "t": [
      0,
      0.1
    ],
    "n": 3,
    "max_step": 0.05
  },
  "stat": {
    "tol": 1e-10
  },
  "chk": {
    "exact": {
      "u": "[x**2*(t + 1)**2,-2*x*y*(t + 1)**2]",
      "p": "(t + 1)**2*(x + y - 1)"
    },
    "tol": 4.0
  },
  "asm": [
    "Manufactured incompressible solution; source generated from the strong PDE",
    "Continuous P2 velocity and P1 pressure; pressure has zero volume mean",
    "Full Dirichlet velocity; nonlinear terms use element quadrature",
    "Use Equilibrium for stationary examples; run for transient examples"
  ]
}
)PHZMODEL"},
        {"pnp", R"PHZMODEL({
  "ver": 1,
  "nm": "Poisson-Nernst-Planck: neutral closed system",
  "dim": 2,
  "fld": ["cp", "cm"],
  "alg": [{"nm": "psi", "a": "eps", "src": "cp-cm"}],
  "par": {"D": 0.1, "eps": 0.5},
  "eq": {"cp": "D*lap(cp)+D*divg(cp,psi)", "cm": "D*lap(cm)-D*divg(cm,psi)"},
  "ini": {"cp": "1+0.05*cos(pi*x)", "cm": "1-0.05*cos(pi*x)"},
  "dom": {"type": "box", "n": [6, 6]},
  "disc": {"type": "fe", "q": 4},
  "run": {"t": [0, 0.5], "n": 11, "method": "BDF", "rtol": 1e-7, "atol": 1e-9},
  "chk": {"mass": ["cp", "cm"], "nonneg": ["cp", "cm"], "tol": 1e-6},
  "asm": ["Dimensionless monovalent electrolyte", "Constant permittivity", "No boundary flux", "Homogeneous Neumann Poisson equation", "Mean-zero potential gauge", "Global charge neutrality is required and checked"]
}
)PHZMODEL"},
        {"rd", R"PHZMODEL({
  "ver": 1,
  "nm": "Reaction-diffusion: OUTPUT reference",
  "dim": 2,
  "fld": ["u", "v"],
  "par": {"Du": 0.1, "Dv": 0.2, "a": 1.0, "b": 0.5, "Lx1": 1.0, "Lx2": 1.0},
  "eq": {"u": "Du*lap(u)+a*u-b*u*v", "v": "Dv*lap(v)+b*u*v-v"},
  "ini": {"u": "2+0.2*cos(pi*x/Lx1)*cos(pi*y/Lx2)", "v": "1+0.1*cos(pi*x/Lx1)"},
  "bc": {"u": {"type": "neumann", "val": 0}, "v": {"type": "neumann", "val": 0}},
  "dom": {"type": "box"},
  "disc": {"type": "sp", "n": [4, 4], "q": 20},
  "run": {"t": [0, 5], "n": 81, "method": "BDF", "rtol": 1e-8, "atol": 1e-10},
  "chk": {"nonneg": ["u", "v"], "tol": 1e-6},
  "asm": ["Dimensionless reaction-diffusion system", "Constant scalar diffusivities", "No boundary flux", "Cosine Galerkin truncation with modes 0..3 on each axis"]
}
)PHZMODEL"},
        {"rkheat", R"PHZMODEL({
  "ver": 1,
  "nm": "Heat modes with exact scheme analysis",
  "dim": 2,
  "fld": ["u"],
  "par": {"D": 0.2, "Lx1": 1, "Lx2": 2},
  "eq": {"u": "D*lap(u)"},
  "ini": {"u": "1+0.2*cos(pi*x)+0.1*cos(2*pi*y/Lx2)"},
  "disc": {"type": "sp", "n": [4, 4]},
  "run": {"method": "RK4", "t": [0, 0.5], "n": 6, "max_step": 0.05},
  "chk": {"mass": ["u"], "exact": {"u": "1+0.2*exp(-D*pi**2*t)*cos(pi*x)+0.1*exp(-D*(2*pi/Lx2)**2*t)*cos(2*pi*y/Lx2)"}, "tol": 1e-6},
  "asm": ["Constant diffusivity and homogeneous Neumann boundaries", "Finite cosine Galerkin space with known exact eigenvalues", "Scheme stability is proved by symbolic algebra in exact arithmetic; no Lean proof or floating-point error certificate is claimed"]
}
)PHZMODEL"},
        {"sphere", R"PHZMODEL({
  "ver": 1,
  "nm": "Surface diffusion on a triangulated sphere",
  "dim": 2,
  "fld": ["u"],
  "par": {"D": 0.1},
  "eq": {"u": "D*lap(u)"},
  "ini": {"u": "1+0.2*z"},
  "dom": {"type": "sphere", "n": 2},
  "disc": {"type": "fe", "q": 4},
  "run": {"t": [0, 0.5], "n": 11, "method": "BDF", "rtol": 1e-8, "atol": 1e-10},
  "chk": {"mass": ["u"], "nonneg": ["u"], "tol": 1e-6},
  "asm": ["Scalar Laplace-Beltrami operator on a fixed polyhedral surface", "The triangulated surface approximates a unit sphere", "No evolving geometry or vector surface mechanics"]
}
)PHZMODEL"},
        {"problem_ns2", R"PHZMODEL({
  "dim": 2,
  "fld": {
    "u": [
      2
    ],
    "p": []
  },
  "par": {
    "nu": "1/100"
  },
  "coef": {
    "f": [
      2
    ]
  },
  "pde": {
    "u": "dt(u)+adv(u,u)=-grad(p)+nu*lap(u)+f",
    "p": "div(u)=0"
  },
  "bc": {
    "u": {
      "type": "dirichlet",
      "val": "[0,0]"
    }
  }
}
)PHZMODEL"},
        {"problem_ns3", R"PHZMODEL({
  "dim": 3,
  "fld": {
    "u": [
      3
    ],
    "p": []
  },
  "par": {
    "nu": "1/100"
  },
  "coef": {
    "f": [
      3
    ]
  },
  "pde": {
    "u": "dt(u)+adv(u,u)=-grad(p)+nu*lap(u)+f",
    "p": "div(u)=0"
  },
  "bc": {
    "u": {
      "type": "dirichlet",
      "val": "[0,0,0]"
    }
  }
}
)PHZMODEL"},
    };
    return d;
}

const std::map<std::string, std::string>& referenceDocs() {
    static const std::map<std::string, std::string> d = {
        {"DSL", R"PHZDOC(# Язык PDE и слабая форма

Версия 0.6.0 ввела ASCII-язык с явными типами полей. Lark разбирает выражения в неизменяемые узлы `CN`; каждый узел хранит оператор, форму тензора и аргументы. Выполнение Python из строки не используется. В версии 0.7.0 добавлен отдельный численный путь `disc.type=th` для смешанных задач скорости и давления; он описан в [FLOW.md](FLOW.md). Ниже ограничения скалярного моста относятся к `disc.type=fe/sp`.

## Ввод

```json
{
  "dim": 2,
  "fld": {"u": [2], "p": []},
  "par": {"nu": "1/100"},
  "coef": {"f": [2]},
  "pde": {
    "u": "dt(u)+adv(u,u)=-grad(p)+nu*lap(u)+f",
    "p": "div(u)=0"
  },
  "bc": {"u": {"type": "dirichlet", "val": "[0,0]"}}
}
```

`[]` означает скаляр, `[2]` — двухкомпонентный вектор, `[2,2]` — матрицу. Список `"fld": ["u", "w"]` задаёт скалярную систему. `par` содержит постоянные численные параметры, `coef` объявляет заданные функции координат и времени без определения их значений. Каждому полю соответствует одно уравнение той же тензорной формы. Мост `fe/sp` вычисляет скалярные поля без свободных функций `coef`. Путь `th` вычисляет скорость и давление, а значения `coef` получает из `src`; см. [FLOW.md](FLOW.md).

Литералы, включая длинные десятичные числа, читаются точно. В символьных постановках параметр можно задавать строкой-дробью. Численные модели сохраняют существующий формат числовых параметров JSON. Допускаются арифметика, `**`, круглые скобки, `exp`, `log`, `sin`, `cos`, `sqrt`, координаты `x,y,z`, время `t`, константа `pi`, вложенные массивы. Отрицательная степень и унарный минус следуют обычному приоритету Python: `-2**2=-4`, `(-2)**2=4`.

| Оператор | Смысл |
|---|---|
| `dt(u)` | Частная производная по времени |
| `grad(u)` | Последний индекс — направление пространственной производной |
| `div(A)` | Свёртка последнего индекса с направлением производной |
| `lap(u)` | Покомпонентный лапласиан |
| `adv(a,u)` | `(a.grad)u`; первый аргумент — пространственный вектор |
| `dot(A,B)` | Последний индекс A сворачивается с первым индексом B |
| `inner(A,B)` | Полная свёртка одинаковых форм; поля вещественные |
| `outer(A,B)` | Внешнее произведение |
| `trans(A)`, `sym(A)` | Транспонирование и `(A+trans(A))/2` |
| `eye()` | Единичный пространственный тензор |
| `divg(a,u)` | Сокращение для `div(a*grad(u))` |
| `intg(u)` | Глобальный интеграл по всей области, сохраняющий тензорную форму |

Обычное `*` требует хотя бы один скалярный аргумент. Формы проверяются до получения слабого выражения. Поддерживаются поля ранга 0–2 и промежуточные выражения ранга до 3. Нет произвольного индексирования, функций Python, атрибутов или доступа к файлам из уравнения. Не поддержаны общие ядра Фредгольма, свёртки и дробные производные.

## Вывод слабой формы

```powershell
.\.venv\Scripts\python.exe -m phz weak problems\ns2.json -o generated\ns2_weak.json
.\.venv\Scripts\python.exe -m phz weak problems\ns3.json -o generated\ns3_weak.json
```

Сохраняются JSON с типизированными деревьями и соседний Markdown с выражениями. Вывод строится из остатка `LHS-RHS=0`, умноженного на тестовую функцию. После интегрирования по частям сумма объёмного и граничного интегралов должна равняться нулю. Для Навье–Стокса с постоянной вязкостью результат эквивалентен:

```text
int_Omega(dt(u).v + adv(u,u).v + nu*inner(grad(u),grad(v))
          - p*div(v) - f.v)
+ int_boundary((p*n - nu*dot(grad(u),n)).v) = 0

int_Omega(q*div(u)) = 0
```

Дополнительно поддерживается запись через напряжения `div(2*nu*sym(grad(u))-p*eye())`. Её граничный член отличается от записи с векторным лапласианом; система выводит член именно из заданного оператора. `div(u)` в уравнении несжимаемости остаётся в объёмном интеграле. Для произведения `a(x)*lap(u)` дифференцируется также коэффициент: объёмный член содержит `grad(a*v)`. Он не подменяется формой оператора `div(a*grad(u))`.

`dirichlet` ограничивает пробное поле и зануляет тестовую функцию на указанной части границы. `neumann` или `flux` задаёт весь поток, спаренный с тестом в остатке. Для `dt(u)=div(D*grad(u))+f` это `q=-D*grad(u).n`. `robin` задаёт `q=coef*(u-val)`. `traction` задаёт противоположный знак, подходящий для внешнего вектора напряжений. Частичные условия используют `on`. Пересечение разных именованных участков проверяется по сетке при расчёте, а символьный вывод предполагает непересекающиеся участки.

Команда `weak` сохраняет незаданные граничные члены явно. Численный запуск применяет существующее соглашение: незаданный внешний поток равен нулю. В его `weak.json` такие участки помечены `default_zero_flux`. Для производных высокого порядка одного условия Дирихле недостаточно: здесь показан один шаг интегрирования по частям, а не автоматически замкнутая постановка произвольного порядка.

## Связь с расчётом

В численной модели вместо `eq` можно указать `pde`. Должна присутствовать одна производная собственного поля с единичным коэффициентом: `dt(u)`. Разрешены локальные реакции, глобальные интегралы и линейные комбинации консервативных операторов `lap(field)`, `div(a*grad(field))` с постоянными внешними множителями. Внутренний коэффициент `a` может нелинейно зависеть от полей. Численный мост не принимает `x*lap(u)`, `lap(u*u)`, вложенные дифференциальные операторы, векторные поля или смешанные алгебраические ограничения. Их отдельный символьный вывод не означает наличия численного решателя.

```powershell
.\.venv\Scripts\python.exe -m phz run examples\dsl2.json
.\.venv\Scripts\python.exe -m phz run examples\int3.json
.\.venv\Scripts\python.exe -m phz study examples\int3.json --levels 2
```

`dsl2` проверяет нелинейную диффузию и одновременно Дирихле, Нейман, Робин. `int3` рассчитывает объёмную модель:

```text
dt(u) = D*lap(u) + b*(intg(u)/intg(1)-u)
```

При нулевом потоке она сохраняет массу. `intg(1)` вычисляется по реальной сетке и даёт объём области. В FEM выражение под интегралом вычисляется в квадратурных точках каждого элемента, затем взвешенные значения суммируются по области. Нелинейные значения узлов не подставляются вместо элементных интегралов. После этого нелокальная величина участвует в обычной сборке нагрузки. Квадратурная чувствительность и чувствительность к сетке остаются отдельными проверками.

При глобальных интегралах стационарный расчёт использует плотный численный якобиан с существующим пределом 256 свободных коэффициентов. Локальная разреженная сборка явно отклоняется: она не учитывает глобальную связь степеней свободы. Экспорт полного полиномиального косинусного `du_dt` для `intg` пока отсутствует; выполняемый расчёт и снимок кода сохраняются в архиве.

Архив содержит исходное `pde` в `model.json`, исходные строки, типизированные деревья, слабую форму, граничные члены, допущения и преобразованные правые части в `weak.json`, а также `weak.md`. Оба новых файла включены в контрольные суммы. Повторение запуска использует сохранённый исходный код и зафиксированные зависимости, включая Lark.

## Основания доверия и пределы

Независимые тесты сравнивают объёмный и граничный члены Навье–Стокса с компонентными формулами в 2D/3D. Для переменных коэффициентов тождество Грина проверяется точным интегрированием полиномов по объёму и всем граням. Сборка `intg(u*u)` проверяется по известным интегралам; дополнительно проверяются масса, аналитическое затухание косинусной моды, изготовленное нелинейное решение и сохранённый replay.

Это символьный вывод и численная верификация, без сертификата Lean/Coq. Начиная с 0.7.0 поддерживаемые смешанные постановки Навье–Стокса рассчитываются через [отдельный P2/P1-путь](FLOW.md); универсальность для произвольных векторных систем этим не установлена. Вспомогательное компонентное раскрытие `dCmp` представляет `intg` интегралом по формальному прямоугольнику с границами `a0,b0,...`; численный `intg` использует действительную область сетки.
)PHZDOC"},
        {"FLOW", R"PHZDOC(# Смешанный расчёт несжимаемого потока

Версия 0.7.0 соединяет типизированную слабую форму с численной сборкой скорости и давления в 2D/3D. Поддерживаются стационарные задачи Стокса и нелинейные задачи Навье–Стокса, а также нестационарный расчёт методом Backward Euler. Геометрия состоит из прямолинейных треугольников или тетраэдров. Скорость представляется непрерывными квадратичными элементами P2, давление — непрерывными линейными P1.

## Запуск

```powershell
.\.venv\Scripts\python.exe -m phz stat examples\nseq2.json -o generated\nseq2.json
.\.venv\Scripts\python.exe -m phz stat examples\nseq3.json -o generated\nseq3.json
.\.venv\Scripts\python.exe -m phz run examples\nst2.json
.\.venv\Scripts\python.exe -m phz run examples\nst3.json
.\.venv\Scripts\python.exe -m phz study examples\nsref2.json --stationary --levels 3
.\.venv\Scripts\python.exe -m phz study examples\nsref3.json --stationary --levels 3
.\.venv\Scripts\python.exe -m phz study examples\nstime.json --kind time --levels 3
```

В настольном окне эти примеры доступны в общем списке. Для `nseq` и `nsref` выбирайте `Equilibrium`, для `nst` — обычный расчёт. Первое поле показывает модуль скорости, второе — давление. Полные компоненты скорости и их порядок сохранены в архиве; график в 2D использует значения в вершинах геометрии и не отображает все квадратичные вариации внутри элемента.

Модель использует `"disc": {"type": "th", "q": 4}`. Поля объявляются как `"fld": {"u": [2], "p": []}` в 2D или с `[3]` для скорости в 3D. Имена полей могут отличаться. `coef` объявляет заданные функции, а `src` задаёт их локальные выражения, например `"coef": {"f": [2]}`, `"src": {"f": "[1,0]"}`. Параметры допускают точные дробные строки, как и символьный DSL.

Поддерживаемая структура несжимаемости и давления:

```text
dt(u) + adv(u,u) = -grad(p) + nu*lap(u) + f
div(u) = 0
```

Для стационарного Стокса можно писать `0=-grad(p)+nu*lap(u)+f`. Поддержана также форма через симметричные напряжения `div(2*nu*sym(grad(u))-p*eye())`. Временная производная скорости должна иметь единичный коэффициент; уравнение давления должно быть `div(u)=0`, а давление в импульсе входить только через `grad(p)` с указанным знаком. Другие допустимые локальные члены собираются из слабой формы, если после интегрирования по частям остаются только значения и первые пространственные производные. В этом вычислительном пути нет глобальных интегралов и вложенных производных высокого порядка.

## Сборка и решение

Узловые базисы построены по барицентрическим координатам. Библиотечный класс `CPn` строит P1–P4 на отрезках, треугольниках и тетраэдрах. Согласование общих рёбер и граней использует топологию и целые барицентрические веса. В пользовательском смешанном решателе сейчас фиксирована пара P2/P1; скалярный FEM-путь по-прежнему использует P1. Полное переключение всех базисов из ТЗ ещё не реализовано.

Типизированная слабая форма раскрывается по компонентам. Коэффициенты при тестовой функции и её градиенте превращаются в выражения для элементной квадратуры. Из этих же выражений символьно получаются производные по значениям, градиентам и временным производным неизвестных. В переносе учитываются обе вариации:

```text
d((u.grad)u)[w] = (w.grad)u + (u.grad)w
```

Разреженные блоки остатка и якобиана собираются интегрированием. Умножение нелинейных значений узлов не заменяет элементную квадратуру. Для P2/P1 требуется `4 <= q <= 20`; отдельное исследование `--kind quad` повышает порядок квадратуры при неизменной сетке. Начиная с 0.8.0 доступен [автоматический контроль](QUADRATURE.md) через `disc.aq`.

Стационарная задача решается демпфированным Ньютоном. В каждом шаге Backward Euler решается вся связанная система скорости и давления. Давление — алгебраическая неизвестная, и его блок матрицы массы равен нулю. Известные граничные значения входят в полный остаток, включая разность по времени; для свободных коэффициентов берётся соответствующая часть якобиана. `stat.tol` ограничивает абсолютную максимальную компоненту слабой невязки и среднее давление. `run.max_step` ограничивает шаг, разбиваемый так, чтобы выходные временные срезы совпадали между исследованиями. `rtol/atol` не управляют точностью данного фиксированного шага; используется `stat.tol`.

## Границы и давление

Пока поддержана скорость Дирихле на всей границе. Можно разделять её на именованные участки, но они должны покрывать всю границу без наложения граней; значения в общих узлах должны совпадать. Скорость может зависеть от времени. При стационарном запуске временная зависимость заданных функций и границ отклоняется.

Давление нормируется условием нулевого среднего. Дополнительная строка и столбец снимают постоянную неопределённость давления. После решения проверяются исходные уравнения несжимаемости, без добавленного множителя: это не позволяет скрыть несовместимый суммарный поток через границу. Для пространств давления до 128 коэффициентов дополнительно проверяется численный ранг матрицы дивергенции на свободных скоростях. Это диагностика лишних нулевых мод, не доказательство inf-sup-оценки. На более крупных сетках эта проверка помечается как невыполненная; невязки линейного и нелинейного решений проверяются всегда.

Начальная скорость нестационарного расчёта проецируется в дискретно несжимаемое пространство при заданных граничных значениях. L2-норма изменения сохраняется. Начальное давление берётся из `ini` с вычитанием среднего; последующие значения находятся совместным решением. Аналитическое давление для сравнения также приводится к нулевому среднему в используемой квадратуре.

## Что проверяется

- Невязки свободных слабых уравнений, граничных значений и среднего давления.
- Моменты `div(u)` по базису давления и отдельно квадратурная L2-норма дивергенции. P2/P1 обеспечивает слабую несжимаемость, но не тождественное обращение дивергенции в ноль в каждой точке.
- L2-ошибки скорости и давления относительно заданных аналитических выражений.
- Кинетическая энергия по квадратуре. Проверка её убывания выполняется только при `chk.energy=true` и не предполагается автоматически для потока с внешней силой или работой границы.
- Сеточная и временная чувствительность по разностям решений; давление везде имеет одну нормировку.

Полиномиальные изготовленные решения проверяют точную воспроизводимость P2/P1, знаки давления, нелинейную сборку и перенос заданных границ во времени. Для более сложных гладких решений отдельно проверяется уменьшение ошибок при сгущении сетки в 2D/3D. Малая слабая невязка на грубой сетке сама по себе не означает точного поля.

Полный контрольный прогон из корня проекта:

```powershell
.\.venv\Scripts\python.exe -m pytest -q -W error --junitxml=generated/tests.xml
.\.venv\Scripts\python.exe -m pytest tests\test_flow.py -q -W error --junitxml=generated/flow_tests.xml
.\.venv\Scripts\python.exe tools\vfy.py
```

Последняя команда проверяет исходные модели вместе со смешанными примерами, запускает исследования уточнения и повторения сохранёнными копиями кода. Итоговые числа записываются в `generated/validation.json` и `generated/VALIDATION.md`. Для 3D-сеточного примера используются три уровня: 402, 2312 и 15468 коэффициентов. Это наиболее затратная часть нового контрольного набора.

## Архив и пределы

Архив содержит исходную постановку, слабую форму/IR, сетку, координаты и связности обоих пространств, порядок компонентов и коэффициентов, свободные и заданные степени свободы, квадратуру, матрицы массы/жёсткости/дивергенции, последний якобиан, состояния, нормы ошибок, энергию, дивергенцию и журнал Ньютона. Названия столбцов журнала записаны в `manifest.json`. Код, зависимости и новые файлы включены в контрольные суммы; обычные `verify`, `replay` и открытие архива работают и для смешанной задачи.

Для начального состояния нестационарного расчёта невязка шага BE ещё не определена. Массив `weak_residual_evaluated` явно отмечает этот срез как непроверенный; нулевой заполнитель в `weak_residual` не участвует в проверке и не является свидетельством решения начального уравнения импульса.

Текущие ограничения: полные границы Дирихле для скорости, фиксированная P2/P1-пара, локальные последовательные вычисления, до 30000 суммарных коэффициентов и дополнительные ограничения памяти квадратуры. Открытые границы, LES/DNS-бенчмарки турбулентности, устойчивость больших чисел Рейнольдса, MPI/GPU, адаптивная сетка, формальное доказательство и гарантированный интервал ошибки этим этапом не установлены. Стационарная сходимость Ньютона не классифицирует динамическую устойчивость.
)PHZDOC"},
        {"GOAL", R"PHZDOC(# Техническое задание на Phyzo

0\. Общие положения

0.1. Назначение системы.

Phyzo — платформа полного цикла для численного исследования PDE: от постановки уравнения до формально верифицированного, воспроизводимого результата, пригодного для публикации и промышленного применения.



0.2. Целевая аудитория.



Research-математики и физики



Инженеры-вычислители (CFD, механика, климат)



ML-исследователи (neural operators, PINN)



Разработчики, встраивающие решатель в свои продукты



0.3. Ключевой принцип.

Система не «помогает считать». Система берёт на себя весь цикл: понимание постановки → вывод схемы → расчёт → доказательство корректности → упаковка результата.



0.4. Метрика успеха проекта (глобальная).

Phyzo решает хотя бы одну задачу, которую не решил ни один существующий открытый инструмент, и результат верифицирован формально.



1\. Символьное ядро (фундамент)

1.1. Ввод PDE на языке, близком к математическому

Требование. Пользователь пишет уравнение в LaTeX-подобном или UFL-подобном DSL:



text

∂u/∂t + u·∇u = -∇p + ν Δu

∇·u = 0

Реализация. Парсер (PLY/ANTLR/Lark) → AST → внутреннее представление (IR) в виде символьного дерева.



Зачем. Без этого всё остальное — ручная работа, а значит, не масштабируется.



Критерий приёмки. Парсятся: линейные и нелинейные PDE, системы, интегро-дифференциальные уравнения, уравнения с дробными производными (опционально).



1.2. Автоматический вывод слабой формы

Требование. Из сильной формы система сама строит слабую форму (умножение на тестовую функцию, интегрирование по частям, учёт граничных условий).



Реализация. Символьные правила + правило Грина + библиотека тождеств (Стокс, Гаусс–Остроградский).



Зачем. Это то, что делает FEniCS/UFL. Без этого нет универсальности.



Критерий приёмки. Для уравнения Навье–Стокса система выдаёт корректную слабую форму без правок пользователя.



1.3. Распознавание симметрий и законов сохранения

Требование. Система находит симметрии Ли и связанные законы сохранения (теорема Нётер).



Реализация. SymPy + реализация алгоритма Бауманна/Ануфриева или подключение внешних CAS.



Зачем. Сохраняющие схемы (symplectic, discrete gradient) без этого не построить, а они критичны для долговременного моделирования.



Критерий приёмки. Для уравнения Кортевега–де Фриза система находит бесконечную иерархию законов сохранения.



1.4. Формальный анализ устойчивости и сходимости

Требование. Система автоматически:



определяет порядок точности схемы,



проверяет L-устойчивость / A-устойчивость,



выводит условие CFL,



оценивает апостериорную ошибку.



Реализация. Символьный анализ Фурье (von Neumann), теория Дальквиста, символьный анализ спектра.



Зачем. Это превращает «численный эксперимент» в «доказательство».



Критерий приёмки. Для явной схемы Эйлера на уравнении теплопроводности система сама выдаёт условие устойчивости.



2\. Библиотека дискретизаций

2.1. Произвольные базисы

Требование. Поддержка:



спектральных (Фурье, Чебышёв, Лежандр, сферические гармоники),



конечных элементов (Lagrange P1–P4, Hermite, Nedelec, Raviart–Thomas),



вейвлетов (Daubechies, spline),



спектральных элементов.



Зачем. Разные задачи требуют разных базисов; универсальность = выживание.



Критерий приёмки. Пользователь переключает базис одной директивой, и результат остаётся корректным.



2.2. Автоматический выбор базиса

Требование. Система анализирует геометрию, граничные условия, гладкость решения и предлагает (или сама выбирает) оптимальный базис.



Реализация. Эвристики + обученная ML-модель на бенчмарках.



Зачем. Снижает порог входа и ускоряет расчёты.



Критерий приёмки. На задаче с разрывом решение выбирает DG-метод, на гладкой — спектральный.



2.3. Квадратуры с автоконтролем

Требование. Гаусс, Гаусс–Лобатто, Монте-Карло, квази-Монте-Карло (Sobol, Halton), адаптивные схемы.



Критерий приёмки. Система сама повышает порядок квадратуры, пока ошибка не упадёт ниже заданного порога.



3\. Численные решатели

3.1. Интеграторы по времени

Требование. Явные (RK4, SSP-RK), неявные (BDF1–5, Radau IIA), IMEX, экспоненциальные (ETD, exponential Runge–Kutta), сохраняющие (symplectic, discrete gradient).



Реализация. Обёртки над SUNDIALS/ARPACK или собственная реализация с адаптивным шагом.



Критерий приёмки. Жёсткая система решается BDF без ручной настройки.



3.2. Adaptive Mesh Refinement (AMR)

Требование. Локальное измельчение сетки на основе апостериорной оценки ошибки. Структурированные и неструктурированные сетки, quadtree/octree.



Реализация. Dual-weighted residual, Zienkiewicz–Zhu, gradient-based indicators.



Зачем. Без AMR невозможно считать реальные задачи (ударные волны, фронты, сингулярности).



Критерий приёмки. Задача с пограничным слоем решается с сеткой в 10 раз меньше, чем uniform, при той же точности.



3.3. Параллелизация

Требование.



MPI для распределённых вычислений,



OpenMP/threads для shared memory,



CUDA/HIP для GPU,



domain decomposition (Schwarz, FETI, multigrid).



Зачем. Без этого размер задачи ограничен памятью одной машины.



Критерий приёмки. Линейная масштабируемость до 1000 ядер на типовой задаче.



4\. Анализ решений

4.1. Поиск стационарных состояний

Требование. Newton–Krylov, pseudo-arclength continuation, deflation.



Критерий приёмки. Находятся множественные равновесия, включая неустойчивые.



4.2. Бифуркационный анализ

Требование. Обнаружение и продолжение:



седло-узловых,



Хопфа,



гомоклинических/гетероклинических орбит,



бифуркаций коразмерности 2+.



Реализация. Интеграция с AUTO-07p или собственная реализация.



Критерий приёмки. Для системы Лоренца система строит полную бифуркационную диаграмму.



4.3. Обратные задачи и UQ

Требование.



Байесовский вывод параметров (MCMC, HMC, VI),



surrogate-модели (Gaussian processes, polynomial chaos),



sensitivity analysis (Sobol indices).



Зачем. Реальные задачи всегда содержат неопределённость.



Критерий приёмки. Оценка параметров PDE с доверительными интервалами на синтетических данных.



5\. Machine Learning слой

5.1. Neural operators

Требование. Интеграция FNO, DeepONet, PINN, transformers для PDE.



Зачем. Ускорение расчётов в 100–1000 раз на повторяющихся задачах.



Критерий приёмки. Surrogate предсказывает решение с относительной ошибкой < 1% на тесте.



5.2. Гибридный решатель

Требование. Система сама решает, где применить ML, а где — классический метод. ML-предсказание верифицируется классическим решателем.



Зачем. ML без верификации — это гадание. Гибрид даёт и скорость, и надёжность.



Критерий приёмки. На сложной задаче гибрид быстрее чистого решателя в 10+ раз при той же точности.



5.3. Обнаружение PDE из данных (SINDy, symbolic regression)

Требование. По временным рядам/полям система находит управляющее уравнение.



Зачем. Это меняет парадигму: не «задали уравнение — решили», а «есть данные — нашли закон».



Критерий приёмки. По данным симуляции Навье–Стокса система восстанавливает уравнение с вязкостью в пределах 5%.



6\. Формальная верификация

6.1. Генерация доказательств

Требование. Каждый результат сопровождается машинно-проверяемым доказательством в Lean 4 / Coq:



корректность дискретизации,



оценка ошибки,



устойчивость.



Зачем. Это единственный способ сделать численный результат научно неоспоримым.



Критерий приёмки. Для линейного уравнения теплопроводности система генерирует Lean-доказательство сходимости, которое компилируется без sorry.



6.2. Сертификаты ошибки

Требование. Interval arithmetic + directed rounding для гарантированных границ ошибки.



Критерий приёмки. Пользователь получает не «≈ 1.234», а «∈ [1.2338, 1.2342] с гарантией».



7\. Продуктовый уровень

7.1. Мультиинтерфейсность

Desktop GUI (текущий Tkinter — оставить как legacy),



Web UI (React + backend),



Jupyter-интеграция,



CLI,



REST/gRPC API,



Python SDK.



7.2. Облако и HPC

Требование. Запуск на кластере из браузера, поддержка SLURM, PBS, Kubernetes, AWS/GCP/Azure.



7.3. Плагинная экосистема

Требование. Сторонние разработчики пишут PDE-модули, базисы, решатели. Версионирование, dependency resolution.



7.4. Воспроизводимость на десятилетия

Требование. Каждый архив — контейнер (Docker/Singularity) + provenance-метаданные (RO-Crate, W3C PROV) + зафиксированные версии зависимостей. Через 20 лет архив запускается.



Критерий приёмки. Архив, созданный сегодня, воспроизводится на эмуляторе через 10 лет.



8\. Научный выход

8.1. Бенчмарки

Требование. Публичный набор задач (от линейных до турбулентности), на которых сравниваются методы.



8.2. Публикации

Требование. Минимум 3–5 статей в топ-журналах (JCP, SINUM, SIAM Review) и конференциях (ICML, NeurIPS, CPP, ITP).



8.3. Решение открытой задачи

Требование. Хотя бы одна задача, не решённая существующими инструментами. Кандидаты:



устойчивость решений 3D Навье–Стокса (частичный прогресс),



обратные задачи реального времени в медицине,



климатические модели с гарантиями,



квантовая гидродинамика.





ВАЖНО:



- **Автоматическое доказательство устойчивости и сходимости** дискретизации. То есть система сама выводит: «эта схема L-устойчива, порядок 2, CFL такой-то» — и, в идеале, формально верифицирует это в Coq/Lean/Isabelle.
- **Символьное обнаружение симметрий и законов сохранения** (теорема Нётер) для заданного PDE — и использование их для построения сохраняющих схем (symplectic, discrete gradient, mimetic).
- **Автоматический поиск бифуркаций и продолжение по параметру** (AUTO, MatCont-класс), включая гомоклинические/гетероклинические орбиты.
- **Обратные задачи и UQ**: Bayesian inference для параметров PDE, MCMC/VI, surrogate-модели (Gaussian processes, neural operators).

Уже требует **знания современной математики на уровне research-математика** + инженерии.









выполняй задачу ВАЖНО ПОМНИ ПРО ПРАВИЛО С ВЕНГЕРСКОЙ НОТАЦИЕЙ! ВСЕ МАКСИМАЛЬНО СОКРАЩЕНО (ВКЛЮЧАЯ КЛАССЫ И ТД))PHZDOC"},
        {"MATH", R"PHZDOC(# Discrete model

## Cosine Galerkin

```text
b_k(x) = prod_j cos(pi*k_j*x_j/L_j)
u_h(x,t) = sum_k U_k(t)*b_k(x)
M_kk = prod_j (L_j if k_j=0 else L_j/2)
lambda_k = pi^2*sum_j k_j^2/L_j^2
lap(b_k) = -lambda_k*b_k
cos(p*x)*cos(q*x) = (cos((p+q)*x)+cos((p-q)*x))/2
```

Products are represented by exact rational weights. Duplicate zero-frequency channels are combined. Intermediate products retain all generated frequencies; projection onto the requested mode set happens after the complete polynomial expression is built. `divg(a,b)` uses the identity below for symbolic export:

```text
div(a*grad(b)) = (lap(a*b)+a*lap(b)-b*lap(a))/2
```

Exact symbolic expressions remain parameterised by domain lengths and physical coefficients. Their numerical evaluation uses NumPy functions generated by SymPy. Nonpolynomial spectral calculations use quadrature and projected intermediate fields; they require independent quadrature and mode convergence studies.

## P1 simplices in an embedding space

```text
d = topological_dimension
g = embedding_dimension
E = [p_1-p_0, ..., p_d-p_0]
H = E^T*E
V = sqrt(det(H))/factorial(d)
grad(N_i) = E*inv(H)*grad_ref(N_i)
M_ij = V*(1+delta_ij)/((d+1)*(d+2))
K_ij = V*dot(grad(N_i),grad(N_j))
u_q = sum_i N_i(x_q)*U_i
R_i = sum_elements sum_q w_q*N_i(x_q)*f(u_q,x_q,t)
M*dU/dt = -D*K*U + R(U,t)
```

The same Gram-matrix construction supports volume simplices and tangential gradients on embedded surface simplices. The latter approximates a surface PDE on the supplied polyhedral geometry. It does not implement evolving-surface transport or vector surface mechanics.

Simplex quadrature is a Duffy transform of tensor Gauss-Legendre rules. `q` is the number of points in each reference direction, not a promised polynomial exactness degree. The Jacobian weights are included. Increasing `q` and comparing integrals is a separate operation from mesh refinement.

## Operators and intermediate fields

```text
P(f) = inv(M)*B^T*W*f_q
L(U) = -inv(M)*K*U
D(a,b) = -inv(M)*sum_j G_j^T*W*diag(a_q)*G_j*P(b)
```

Here `B` evaluates basis functions at quadrature points and `G_j` evaluates their derivatives. For an existing discrete field, `P(b)` means its coefficients without reprojecting. For a nonlinear expression it means the L2 projection of that expression. Outer sums and constant factors are assembled linearly. Nested operators therefore specify a particular mixed, projected discretisation. Products of local reaction expressions are evaluated together at quadrature points, without replacing them by products of nodal interpolants.

Each weak differential operator starts with zero natural boundary load. For conservative second-order equations, prescribed total outward fluxes are then assembled once per equation. Nested differential operators retain homogeneous natural data; separate heterogeneous conditions for intermediate fields are not implemented.

## Boundary assembly and time-dependent lifting

```text
du/dt = div(J) + f
q = -J.n
q_N = val
q_R = coef*(u-val)
B_i = -int_boundary(N_i*q)
U_c = g(t)
M_ff*dU_f/dt = R_f(U_f,g,t) + B_f(U_f,g,t) - M_fc*dg/dt
```

`q` includes the complete physical transport coefficient. The implementation does not multiply a prescribed flux by the diffusion constant a second time. Robin coefficients and external values are evaluated at boundary quadrature points; P1 traces supply the field values there. The boundary rule uses the same per-direction Gauss order as the volume rule, with a Duffy transform on triangular faces. In 1D, endpoint measure is one.

Only free coefficients enter the time integrator. Boundary values are reconstructed from their expressions at every RHS evaluation and saved time. Their derivatives are obtained symbolically. Backward Euler is applied to the reduced free ODE, including the exact `dg/dt` at the new time. Its Newton residual is the reduced ODE step residual, extended by zero at constrained nodes for the spatial L2 norm. Boundary forcing can inject mass or energy, so closed-system conservation or decay must not be assumed for arbitrary prescribed data.

Coordinate tags select entire exterior facets at coordinate extrema. Custom tags explicitly list exterior facet vertices. Refinement splits these same facets and preserves their names; replay uses the stored mesh and tags. Overlapping facets within a field are rejected. Shared Dirichlet corner nodes must agree as functions of time after coordinate and parameter substitution.

## Variational model

```text
F = int(f(u_1,...,u_n) + sum_i kap_i*|grad(u_i)|^2/2)
mu_i = df/du_i - kap_i*lap(u_i)
AC: du_i/dt = -mob_i*mu_i
CH: du_i/dt = mob_i*lap(mu_i)
```

`kap_i` and `mob_i` are constant parameters. The local derivative is obtained symbolically. For FEM, the chemical potential is an L2-projected P1 field. The semidiscrete energy and its derivatives use the same quadrature. A generic BDF/Radau time integrator does not imply unconditional energy stability; energy is measured and reported.

## Poisson constraint and electrolyte example

```text
-a*lap(psi) = rho
a*K*Psi = M*P(rho)
int(psi) = 0
int(rho) = 0
d(cp)/dt = D*lap(cp)+D*div(cp*grad(psi))
d(cm)/dt = D*lap(cm)-D*div(cm*grad(psi))
rho = cp-cm
```

The Neumann compatibility condition is checked before each Poisson solve. A mean constraint fixes the null mode. The elliptic solve is performed for every evaluation of the dynamic RHS; the algebraic variables and their residuals are saved. This is elimination of a specific linear elliptic constraint, not a generic high-index DAE solver.

Compatibility uses the physical integral `aMass dot rho_coeff`, not the unweighted sum of modal load entries. The latter is inappropriate for a cosine basis because nonconstant modes have zero integral.

`aMass=M*aOne`, where `aOne` represents the constant unit field: all ones for nodal P1, and only the zero mode for the cosine basis. These integrals are analytic and independent of the selected quadrature order. Projecting a nonlinear or coordinate-dependent source still involves quadrature and must be checked separately.

## Numerical Jacobian

The numeric path converts state vectors to floating point before applying perturbations or prescribed boundary values. Central differences use the actual representable coordinate separation. When a perturbation violates a real-domain constraint, its size is reduced; a valid one-sided difference is retained as a fallback if central perturbations remain unavailable. A non-finite or unresolvable derivative raises an error. The exact polynomial spectral path retains its symbolic Jacobian. Domain-aware differencing does not make a singular derivative finite and does not provide a positivity theorem for the integrator.

The PNP example is nondimensional, uses constant permittivity and unit valences, and has no fixed charge. Other source expressions can include fixed charge if global neutrality is satisfied. The implemented consistent-mass Galerkin/BDF scheme does not inherit positivity or unconditional energy theorems from specialised schemes in the literature.

## Stationary states and constrained linearization

```text
F(U_f) = RHS(U_f,g)_f
A*U_f = b
A*Z = 0
Z^T*Z = I
U_f = U_0 + Z*w
r(w) = Z^T*F(U_0+Z*w)
J_r = Z^T*DF*Z
J_r*dw = -r
w_next = w + alpha*dw
```

Dirichlet values are fixed and time independent. Each row of `A` is the free mass vector of one requested field, normalized to unit Euclidean norm. The prescribed contribution from constrained nodes is subtracted from `b`. Since constraints refer to distinct fields, these rows have disjoint support. The dense solver computes an orthonormal null-space basis by SVD. A compatible starting point is obtained by orthogonal correction of the initial free coefficients.

Newton uses backtracking until the norm of the projected residual decreases sufficiently. Domain-invalid trial values are rejected. Termination requires the infinity norm of the **full free RHS**, before tangent projection, to meet `stat.tol`. Requested mass errors, normalized by the corresponding mass-vector norms, are checked separately. Singular Jacobians and incompatible constraints raise errors rather than producing a successful equilibrium report.

`DF*Z` is differentiated directly along admissible perturbations. For polynomial spectral problems the exact symbolic Jacobian is used. This is essential for Neumann Poisson coupling: arbitrary individual ion-coefficient perturbations would violate neutrality. The stationary coupled implementation therefore requires an affine source with constant field coefficients and an explicit fixed mass for every source field.

The reported eigenvalues belong to `J_r`, which linearizes the semidiscrete **dynamic RHS**; no additional mass-matrix multiplication is applied. Dirichlet perturbations vanish. With mass constraints, the code also reports `norm(A*DF*Z)/max(1,norm(DF*Z))`; excessive leakage prevents a stability classification on that subspace. Stability is local and pertains to the selected finite-dimensional linearization. A single stationary result does not test temporal conservation or energy decay.

The archive stores one equilibrium state at reference `t=0`, the accepted Newton iteration residuals and damping factors, and mass targets. A dense tangent basis is stored only up to 256 free coefficients; the reduced Jacobian is stored when eigenanalysis is performed. Stationary refinement studies keep the first-level mass targets fixed.

## Sparse stationary FEM

For a local reaction `f(u)` and a conservative term `div(a(u,x)*grad(u_p))`, the same volume quadrature gives the residual and its derivative with respect to field `u_j`:

```text
R_i = sum_q w_q*(B_qi*f_q - sum_d G_dqi*a_q*(G_d*U_p)_q)
J_ij = B^T*diag(w*df/du_j)*B
       - sum_d G_d^T*diag(w*da/du_j*(G_d*U_p))*B
       - delta_jp*sum_d G_d^T*diag(w*a)*G_d
J_robin = -B_boundary^T*diag(w_boundary*k)*B_boundary
```

Indices `i,j` in the second expression denote field blocks; `B` is the P1 evaluation matrix and `G_d` its physical or tangential derivative. Local derivatives are symbolic, but their values and products are evaluated at quadrature points. This includes cross-field diffusivity derivatives and does not substitute nodal nonlinearities. The weak Jacobian is sparse before applying the inverse consistent mass matrix. Prescribed boundary coefficients are eliminated from both its rows and columns. Neumann data is independent of the unknown fields and has zero derivative.

The supported sparse form is a sum of local reactions, constant multiples of `lap(field)` and `divg(local,field)`, with no algebraic fields. Nonlocal projected intermediate operators retain the existing dense solver. The sparse assembly never silently changes their discretization.

```text
M_ff*dU_f/dt = R(U_f,g)
J_w = dR/dU_f
[J_w A^T; A 0]*[dU; lambda] = [-R; b-A*U_f]
U_next = U_f + alpha*dU
```

The bordered system is an auxiliary Newton correction system. Its multipliers are not additional physical PDE terms. Constraint correction removes numerical mass drift. Backtracking uses the full weak-residual Euclidean norm; success requires the infinity norm of `M_ff^-1*R` below `stat.tol` and a second check through the existing full RHS evaluator. Thus an augmented solution with nonzero original PDE residual is not accepted. Each sparse LU solve is checked by its relative linear residual; very small LU pivots cause a singularity/rescaling error. This pivot check is a numerical guard, not a condition-number estimate. Cancellation and the time limit are checked between assembly and factorization calls; a running native factorization cannot be interrupted mid-call.

For small sparse solves, stability uses `Z^T*M_ff^-1*J_w*Z`, including the same constraint-invariance check as the dense path. Above 256 free coefficients the full spectrum is explicitly `not_run_size_limit`; no sparse partial-spectrum stability claim is made. Large constrained solves avoid storing an `N x N` dense null-space basis. The protective sparse solve limit is 20,000 free coefficients, additionally bounded by mesh and quadrature limits; memory and time remain problem dependent.

`stationary_weak_jacobian.npz` stores the final free weak Jacobian in a SciPy sparse format. Its checksum, copied assembly code, residual history, mass rows, solver selection and any dense-fallback reason are part of the run archive. The mesh study and quadrature study remain independent of the derivative verification.

## Verification and refinement

The reference output is checked by symbolic subtraction of every generated RHS from every supplied RHS. Numerical quadrature independently checks the spectral projection. Element reaction loads are checked against exact simplex monomial integrals, including a case that distinguishes quadrature assembly from nodal substitution.

Manufactured heat solutions measure absolute L2 error. Nested mesh studies compare transferred coarse solutions and fine solutions at identical saved times using fine-grid quadrature. A finer-grid difference is not labelled as a rigorous error estimator. Geometry, quadrature and time approximation may require their own independent studies.

Time studies use fixed-step backward Euler or explicit FE/RK2/SSPRK3/RK4 with successive step halving and fixed spatial discretisation, quadrature, method, saved times and solver tolerances. Every level stores actual step lengths and all inter-level L2 differences. For at least three levels, `log2(norm(U_h-U_h/2)/norm(U_h/2-U_h/4))` is reported as an empirical order; values below the numerical comparison threshold are omitted. Spatial or quadrature error may remain even when this difference is small. Adaptive `solve_ivp` methods are deliberately excluded from fixed-step studies because their `rtol`/`atol` are local error controls rather than prescribed step sizes.

## Exact symbolic RK analysis

```text
c = A*1
R(z) = 1 + z*b^T*(I-z*A)^(-1)*1 = P(z)/Q(z)
gamma(t) = |t| * product(gamma(child))
Phi(leaf) = 1
Phi(t) = product(A*Phi(child))
b^T*Phi(t) = 1/gamma(t)
```

Unordered rooted trees are generated recursively by nondecreasing forests. The tree counts through order 8 are `1,1,2,4,9,20,48,115`. The elementary-weight conditions, not just agreement of `R(z)` with `exp(z)`, determine classical nonlinear order. The first failed order gives an exact order; passing every checked order gives only a lower bound. All tableau values and conditions use exact rationals. Sufficient smoothness and the usual local-to-global ODE assumptions remain necessary; the order is not a PDE convergence theorem.

For `z=-x+I*y`, `x>=0`, the implementation expands `|Q(z)|^2`, `|det(I-z*A)|^2` and `|Q(z)|^2-|P(z)|^2`. Nonnegative rational coefficients and even powers of `y` certify nonnegativity on the closed left half-plane. A positive constant coefficient in the first two polynomials certifies nonvanishing. The determinant check prevents cancelling a pole of the internal stage equations from incorrectly certifying a usable method. These are sufficient certificates; failure to find one is not a disproof. An exact pole or amplification exceeding one at a negative rational point is a disproof. Otherwise the result is unknown. L-stability requires A-stability and a zero rational-function limit at infinity.

The negative-real interval containing the origin is obtained by exact univariate polynomial inequalities and exclusion of denominator and stage poles. Its endpoint is kept as an exact algebraic expression; the decimal endpoint is explicitly an approximation.

For constant-coefficient scalar heat flow in the finite cosine space, orthogonality and known Laplace eigenvalues give `rho=D*pi^2*sum((n_i-1)^2/L_i^2)`. The RK condition is `dt*rho<=r`, with endpoint inclusion inherited from the interval. This is a semidiscrete L2 perturbation-stability statement in exact arithmetic, including perturbations in modes absent from the initial state. It is not a bound on spatial truncation, initial projection, time error, or rounding. For P1 FEM only the conditional statement in terms of the largest generalized eigenvalue of `(K_ff,M_ff)` is emitted; no floating-point eigensolve is presented as a verified enclosure.

These certificates are checked by symbolic algebra in Python/SymPy. They have not yet been translated into a proof-assistant kernel. The archive states `formal_proof_assistant=false`. Explicit stepping uses the same tableau definitions; time-dependent Dirichlet values are lifted at every stage. Intermediate stored-step interpolation is piecewise linear and must not be confused with an order-p RK dense output.
)PHZDOC"},
        {"QUADRATURE", R"PHZDOC(# Автоматический контроль квадратуры

Версия 0.8.0 добавляет подбор порядка квадратуры для скалярного FEM, косинусного Галёркина и смешанного P2/P1-решателя. Он доступен при обычном расчёте и поиске равновесия, включая настольное окно. Геометрия, базис, исходные функции, моменты вывода и настройки решателя остаются заданными моделью.

## Как включить

Добавьте `aq` внутрь `disc`:

```json
"disc": {
  "type": "fe",
  "q": 3,
  "aq": {
    "atol": 1e-8,
    "rtol": 1e-7,
    "qmax": 10,
    "step": 2,
    "confirm": 2,
    "timeout": 300
  }
}
```

`q` — начальный порядок, `qmax` — верхний предел, `step` — приращение. Последний шаг сокращается, если иначе превышен предел. `confirm` задаёт число последовательных успешных сравнений: от 2 до 4. Двух сравнений требуют как минимум три расчёта. `timeout` ограничивает всю процедуру; ограничения времени отдельного решателя также действуют.

Пустой объект `"aq": {}` включает значения по умолчанию: `atol=1e-8`, `rtol=1e-5`, `confirm=2`, `timeout=300`. Для FEM/P2–P1 предел по умолчанию равен 10 и шаг 2; `qmax` можно увеличить до 20. Для косинусного базиса предел равен 100 и шаг 4. Начальный порядок должен оставлять место для требуемых сравнений. Без `aq` используется обычный расчёт с фиксированной квадратурой.

Для скалярного P1 минимальный порядок равен 2 на отрезках и треугольниках, 3 — на тетраэдрах. Последнее ограничение необходимо для точного интегрирования матрицы массы после преобразования Даффи. Для смешанного P2/P1 минимум равен 4. Ограничения числа квадратурных точек и памяти сохраняются.

## Критерий остановки

На каждом новом порядке выполняется полный расчёт. Для каждой компоненты сравнения используется условие:

```text
err <= atol + rtol * max(norm_low, norm_high)
```

Проверяются:

- L2-разность полей между соседними порядками во все сохраняемые моменты.
- Действие обоих пространственных операторов на **одних и тех же** состояниях обоих расчётов. Отдельно проверяются начальные состояния и начальные приближения Ньютона.
- Алгебраические поля Пуассона, если они есть.
- Энергия на одинаковых состояниях для вариационной модели.

Для скалярного пути действие оператора измеряется L2-нормой полудискретной правой части. Для смешанного потока сравниваются свободные слабые нагрузки в двойственной норме согласованной матрицы массы; скорость и давление контролируются отдельно. В этом сравнении смешанный оператор пространственный: временная производная не подставляется. Сам расчёт во времени и сравнение его полей выполняются полностью.

Абсолютный допуск применяется к перечисленным величинам в единицах модели. Значения `absolute`, `scale` и `ratio` в отчёте — максимумы по проверенным состояниям; эти максимумы могут достигаться в разные моменты. Для каждой отдельной оценки отношение вычисляется до взятия максимума.

`stat.mass` с числовыми значениями сохраняется. Значение `initial` каждый раз вычисляется по той же исходной функции и текущей квадратуре; соответствующие массы записываются в историю. Это позволяет уточнить неточно вычисленную начальную массу, а не закреплять ошибку первого порядка. При отдельном исследовании сетки действует прежняя политика фиксации выбранной массы.

## Результат и воспроизведение

Если все сравнения прошли нужное число раз подряд, проверка получает `pass`. Если достигнут `qmax`, а подтверждений недостаточно, последний расчёт сохраняется с `quadrature_control: fail`. Даже успешная последняя пара не отменяет требование последовательных подтверждений. CLI возвращает код 2 для такого результата; отмена, превышение времени и ошибки решателя прерывают процедуру явно.

Архив содержит:

- `model.json` — фактическую модель с выбранным фиксированным порядком.
- `quad_model.json` — запрос на полный автоматический подбор и сохранённую геометрию.
- `quad.json` — допуски, порядки, историю норм, статусы и контрольные суммы промежуточных состояний.
- Обычные поля, матрицы, проверки, код и зависимости.

`replay.py` и команда `replay` повторяют всю процедуру подбора. При сравнении чисел следует учитывать округление: контрольная сумма идентифицирует конкретный массив, но не обещает побитового совпадения между процессами и платформами.

Для независимого исследования сетки или шага времени используйте сохранённый `model.json`. Одновременное изменение квадратуры исказило бы смысл такого исследования, поэтому `study` явно отклоняет модель с активным `aq`.

## Готовые примеры

В настольном окне откройте пример и нажмите `Equilibrium`:

| Файл | Постановка |
|---|---|
| `aq1.json` | Экспоненциальная нагрузка, P1 в 1D |
| `aq3.json` | Та же нагрузка, P1 в 3D |
| `aqsp.json` | Косинусное пространство в 1D |
| `aqns2.json` | Нелинейный несжимаемый поток в 2D |
| `aqns3.json` | Нелинейный несжимаемый поток в 3D |

В примерах потока порядок увеличивается на 1, а для 3D задан `qmax=12`: прежнего предела 10 недостаточно для двух последовательных подтверждений при сохранении заданного допуска.

Для `aq1.json` также доступен обычный расчёт во времени.

```powershell
.\.venv\Scripts\python.exe -m phz stat examples\aq3.json -o generated\aq3.json
.\.venv\Scripts\python.exe -m phz run examples\aq1.json
.\.venv\Scripts\python.exe -m pytest tests\test_aq.py -q -W error --junitxml=generated/aq_tests.xml
.\.venv\Scripts\python.exe tools\vqa.py
```

## Смысл проверки

Это эмпирический контроль чувствительности к порядку квадратуры. Он не является гарантированной оценкой ошибки PDE, не проверяет все промежуточные моменты времени и не заменяет сгущение сетки. Узкий пик или осцилляции, пропущенные несколькими правилами, могут дать ложное согласие. Требование последовательных подтверждений уменьшает отдельные риски, но не устраняет их математически. Ограничения оценивания рассмотрены в [обзоре P. Gonnet](https://arxiv.org/abs/1003.4629).

Точность решателя должна позволять различать запрошенные изменения. Слишком строгий допуск квадратуры при грубом решении нелинейной системы может привести к недостижению критерия. Правила Гаусса–Лобатто, MC/QMC, локальное разбиение интеграла и гарантированные интервальные оценки пока не добавлены.
)PHZDOC"},
        {"RANDOM_TESTS", R"PHZDOC(# Протокол 1000 случайных проверок

Дата: 28 сентября 2026. Базовый seed: `20260928`.
Протокол фиксируется до первого запуска. Исходный вычислительный код не меняется:
после испытаний пользователь отдельно решает, какие дефекты исправлять.

Это 1000 отдельных случаев pytest, а не 1000 утверждений внутри одного теста.
Каждый случай имеет независимый поток NumPy PCG64, полученный через
`SeedSequence([seed, номер_семейства, номер_случая])`.
Параметры генерируются заново из этого ключа; случай не отбрасывается при неудаче.
Внутри семейств перемешан порядок номеров. Seed можно переопределить через
`PHZ_SEED`; другой seed означает другую кампанию, а не повтор исходной.

| Семейство | Случаев | Независимое основание проверки |
|---|---:|---|
| Геометрия | 120 | Объём через сингулярные числа, точные моменты симплекса, энергия линейного поля |
| Базисы P1–P4 | 120 | Воспроизведение полинома и его градиента, точный интеграл квадрата |
| Нелинейная FEM-сборка | 140 | Аналитические интегралы степеней барицентрических координат для нагрузки и Якобиана |
| Спектральная проекция | 100 | Независимое дискретное косинусное интегрирование на равномерных средних точках |
| Эволюция | 100 | Затухание тепловых мод и матричная экспонента для связанной линейной системы |
| Границы | 80 | Смешанные условия для аффинного точного решения; баланс потока и массы |
| Равновесия | 60 | Заранее заданный корень нелинейной реакции и постоянное решение с заданной массой |
| Смешанный поток | 40 | Точное течение Пуазейля P2/P1 в 2D/3D, давление, дивергенция |
| Сгущение | 60 | Сохранение интегралов при переносе и убывание ошибки синусоидальной задачи Пуассона |
| Разбор выражений | 80 | Независимая арифметика Fraction и численные функции math |
| Некорректные входы | 60 | Отказ при недопустимой геометрии, выражениях, параметрах запуска |
| Архивы | 40 | Сохранение модели и полей, повторный запуск части архивов, обнаружение повреждений |
| **Всего** | **1000** | |

Размерность, ориентация, форма и масштаб симплексов, длины областей,
коэффициенты, знаки и амплитуды полей меняются независимо. Для геометрических
проверок масштаб рёбер лежит в пределах 0.01–100, отношение заданных главных
масштабов — не более 4. Отдельно проверяются вложенные сегменты и треугольники.
Число элементов и шаг времени ограничены, чтобы набор был пригоден для
регулярного запуска на рабочем компьютере.

Эталон нелинейной сборки не вызывает квадратуру или сборщик проекта:
`int(prod(lambda_j**a_j)) = d!*V*prod(a_j!)/(d+sum(a_j))!`.
Эталон спектрального произведения не использует символьное умножение мод
проекта: число равномерных точек достаточно для точной дискретной
ортогональности всех частот произведения и тестовой моды.

Допуски заданы до запуска: геометрия и полиномиальные интегралы — относительные
погрешности порядка 1e-10–1e-8 с масштабом объёма/градиента; спектральная
проекция — 2e-10; точные стационарные решения — 2e-7; смешанный поток — 2e-7
для скорости, 2e-6 для давления. Для временных схем сравнение идёт с
непрерывным аналитическим решением; допуск зависит от порядка метода и шага,
для адаптивных методов — 3e-7. При сгущении вдвое ожидается уменьшение L2
ошибки гладкой задачи Пуассона как минимум в 2.5 раза, для P1 теоретический
асимптотический порядок — 2. Все точные значения допусков содержатся в тестах.

Если тест падает, исходный случай, seed и первый результат сохраняются.
Исправлять ошибки самого тестового стенда допустимо только с отдельным
объяснением; ослаблять математическое требование ради прохождения нельзя.
Успешный прогон не доказывает корректность произвольной PDE, физическую
адекватность модели, качество GUI или работоспособность на больших сетках.

Запуск из корня проекта:

```powershell
.\.venv\Scripts\python.exe -m pytest tests/test_random.py -q -W error -o junit_family=legacy --junitxml=generated/random_tests.xml
```

Один случай повторяется по его полному идентификатору из отчёта, например:

```powershell
.\.venv\Scripts\python.exe -m pytest "tests/test_random.py::test_geo[000]" -q -W error
```
)PHZDOC"},
        {"SOURCES", R"PHZDOC(# Исследованные источники

Источники использованы для выбора математической структуры и проверок. Они не являются доказательством корректности конкретного кода; такое основание дают независимые тесты реализации и исследования сходимости. Полные тексты статей в проект не копировались.

1. E. Tadmor, **The Spectral Fourier Approximation** и **Aliasing and spectral accuracy**. Проекция и отличие спектральной аппроксимации от интерполяционной обработки нелинейностей. В реализации полиномиальный косинусный случай вычисляется точными тождествами произведений; неполиномиальный случай отмечается как квадратурный.
   - https://math.umd.edu/~tadmor/lecture-notes/hyperbolic/node28.html
   - https://www.math.umd.edu/~tadmor/lecture-notes/hyperbolic/node16.html

2. A. Logg, K.-A. Mardal, G. N. Wells (eds.), **Automated Solution of Differential Equations by the Finite Element Method**. Автоматизация вариационных постановок, сборка и разделение уровней модели.
   - https://fenicsproject.org/pub/book/book/fenics-book-2011-06-14.pdf

3. H. P. Langtangen, A. Logg, **Solving PDEs in Python: The FEniCS Tutorial I**. Изготовленные решения и проверка зависимости ошибки от размера сетки.
   - https://pub.fenicsproject.org/tutorial/pdf/fenics-tutorial-vol1.pdf

4. G. Dziuk, C. M. Elliott, **Finite element methods for surface PDEs**, Acta Numerica 22 (2013), 289–396. Скалярные поверхностные PDE на триангуляциях. Наша реализация ограничена неподвижной кусочно-плоской геометрией.
   - https://doi.org/10.1017/S0962492913000056

5. **DOLFINx: Cahn–Hilliard equation**. Смешанная постановка через концентрацию и химический потенциал; вариационная производная локальной свободной энергии.
   - https://docs.fenicsproject.org/dolfinx/main/python/demos/demo_cahn-hilliard.html

6. M. L. Huggins, **Thermodynamic Properties of Solutions of Long-Chain Compounds** (1942). Физический контекст решёточной энергии смешения. В `fh.json` реализована заданная в IDEA безразмерная плотность смешения; эластичность сетки не добавляется без отдельной механической постановки.
   - https://doi.org/10.1111/j.1749-6632.1942.tb47940.x

7. **High-order space-time finite element methods for the Poisson-Nernst-Planck equations: Positivity and unconditional energy stability** (preprint 2021; journal 2022). Показывает, что гарантии положительности и энергетической устойчивости требуют специально построенной схемы. Наш пример PNP не заявляет эти гарантии: концентрации и массы проверяются численно.
   - https://arxiv.org/abs/2105.01163
   - https://www.sciencedirect.com/science/article/pii/S0045782522002559

8. **SciPy solve_ivp**. Публичные интерфейсы BDF/Radau, допуски и диагностика интегратора. Внутренний журнал невязок Ньютона через этот интерфейс не предоставляется; собственный BE добавлен для явной регистрации итерационных невязок.
   - https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html

9. **SymPy calculus and numerical computation**. Символьные производные и переход к вычислимым выражениям. Пользовательский ввод предварительно разбирается ограниченным AST-парсером.
   - https://docs.sympy.org/latest/modules/calculus/index.html
   - https://docs.sympy.org/latest/modules/numeric-computation.html

10. **UFL form language**. Поля, вариационные формы и различение размерности области и пространства вложения. UFL не является обязательной зависимостью этой реализации.
    - https://docs.fenicsproject.org/ufl/main/manual/form_language.html

11. **W3C PROV Overview**. Концепция происхождения результата через связи между входами, действиями и производными объектами. Архив Phyzo хранит такие связи в собственном JSON-формате и не заявляет соответствия сериализации PROV.
    - https://www.w3.org/TR/prov-overview/

12. H. P. Langtangen, A. Logg; реализация J. S. Dokken, **Setting multiple Dirichlet, Neumann, and Robin conditions**, FEniCSx tutorial. Слабая форма смешанных условий и интегрирование по помеченным участкам. Соглашение `q=-k*du/dn` использовано для положительного потока наружу; его вклад в динамическую правую часть имеет отрицательный знак. Наша реализация использует собственную сборку P1, а не код решателя DOLFINx.
    - https://jsdokken.com/dolfinx-tutorial/chapter3/robin_neumann_dirichlet.html

13. **Python tkinter: Threading model**. Проверка жизненного цикла окна и доступа к Tk из фоновых вычислений. Рабочий поток передаёт данные через очередь и не удерживает объект интерфейса; освобождение ресурсов окна выполняется в основном потоке.
    - https://docs.python.org/3/library/tkinter.html#threading-model

14. **SciPy differentiate.jacobian**. Контекст численного дифференцирования: размер шага и направление разности являются частью вычислительного метода. В Phyzo используется собственная ограниченная процедура уменьшения возмущения при нарушении области определения; вызова scipy.differentiate.jacobian в ядре нет.
    - https://docs.scipy.org/doc/scipy/reference/generated/scipy.differentiate.jacobian.html

15. **PETSc SNES: Nonlinear Solvers**. Структура нелинейной задачи, метод Ньютона и поиск длины шага по невязке. В Phyzo реализованы собственные плотный и разреженный пути Ньютона; PETSc не является зависимостью.
    - https://petsc.org/release/manual/snes/

16. **SciPy linalg.null_space**. Ортонормированный базис нуль-пространства через SVD. Используется для параметризации возмущений, сохраняющих заданные массы при стационарном решении.
    - https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.null_space.html

17. **FEniCS: Nonlinear Poisson equation**. Линеаризация слабой формы диффузии с коэффициентом `1+u**2` содержит как производную градиента, так и производную коэффициента. Использовано при выводе собственного якобиана, сверенного численным дифференцированием и точными интегралами; DOLFIN не является зависимостью.
    - https://olddocs.fenicsproject.org/dolfin/2019.1.0/python/demos/nonlinear-poisson/demo_nonlinear-poisson.py.html

18. **SciPy 1.15.3: splu**. Разреженное LU-разложение через SuperLU, формат CSC и решение факторизованной системы. Применяется к матрице Ньютона с дополнительными строками и столбцами ограничений массы.
    - https://docs.scipy.org/doc/scipy-1.15.3/reference/generated/scipy.sparse.linalg.splu.html

19. F. Bornemann, **Runge-Kutta methods, trees, and Mathematica** (2002). Рекурсивные корневые деревья и полные условия порядка Бутчера. В Phyzo реализована собственная рекурсия на Python; совпадение ряда функции устойчивости отдельно не считается доказательством нелинейного порядка.
    - https://arxiv.org/abs/math/0211049
    - https://arxiv.org/html/math/0211049v1

20. **SciML RootedTrees.jl**. Независимый открытый проект с формулировкой элементарных весов и плотности дерева `gamma`; используется как источник определения, не как зависимость Python-ядра.
    - https://github.com/SciML/RootedTrees.jl

21. MIT, **18.336 Numerical Methods for Partial Differential Equations**. Связь анализа фон Неймана со спектром оператора шага для линейной однородной задачи. В текущей реализации ограничение выводится для фактического конечного косинусного пространства; формула равномерной конечно-разностной сетки не подставляется вместо FEM.
    - https://math.mit.edu/~stevenj/18.336/

22. **UFL 2025.2: Form language**. Тензорные формы, соглашения `grad`, `div`, `dot`, `inner` и разделение коэффициентов и тестовых функций. В Phyzo используются собственные Lark-грамматика и типизированное дерево; UFL не является зависимостью.
    - https://docs.fenicsproject.org/ufl/2025.2.0.post0/manual/form_language.html

23. J. S. Dokken, **The Navier–Stokes equations**. Объёмные и граничные члены слабой формы, соглашение о градиенте скорости и напряжениях. Использовано для независимого сопоставления знаков давления, вязкости и внешнего потока. Алгоритм расщепления из учебника здесь не реализован.
    - https://jsdokken.com/dolfinx-tutorial/chapter2/navierstokes.html

24. **Lark: Grammar reference**. EBNF-грамматики, терминалы, правила и приоритеты; используется LALR-парсер. Зависимость зафиксирована как `lark==1.3.1`.
    - https://lark-parser.readthedocs.io/en/stable/grammar.html

25. **DefElement: Lagrange** и **Degree 2 Lagrange on a tetrahedron**. Полиномиальные пространства и узловые степени свободы на симплексах; использованы для сверки собственной барицентрической реализации P1–P4. Узловое свойство, воспроизведение полиномов, градиенты и интегралы проверяются независимо. Basix не является зависимостью.
    - https://defelement.org/elements/lagrange.html
    - https://defelement.org/elements/examples/tetrahedron-lagrange-equispaced-2.html

26. **DOLFINx: Stokes equations using Taylor-Hood elements**. Смешанные пространства скорости/давления, структура блоков и постоянное нуль-пространство давления. В Phyzo используется собственная P2/P1-сборка и ограничение нулевого среднего; PETSc/MPI из примера здесь не подключены. Знаки проверяются по принятому сильному уравнению Phyzo и изготовленным решениям.
    - https://docs.fenicsproject.org/dolfinx/main/python/demos/demo_stokes.html

27. P. Gonnet, **A Review of Error Estimation in Adaptive Quadrature**. Обзор оценок ошибки и случаев ложного согласия правил. Использован для определения границ утверждений автоконтроля Phyzo: последовательные сравнения дают численное свидетельство, но не гарантированный интервал ошибки. Реализованный глобальный подбор порядка не выдаётся за алгоритм адаптивного разбиения из статьи.
    - https://arxiv.org/abs/1003.4629

28. **SciPy: quad**. Документированный смысл абсолютного и относительного допусков и оценочного характера возвращаемой ошибки. Phyzo использует собственный критерий `atol + rtol*scale` для полей и операторов; QUADPACK не используется для элементной сборки.
    - https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
)PHZDOC"},
        {"STATUS", R"PHZDOC(# Состояние полного ТЗ

Актуальная редакция цели из вложения `ee4bf789-0cb2-456d-b4f2-e7c983136fab` сохранена без сокращений в [GOAL.md](GOAL.md), включая усиленное требование к сокращённой венгерской нотации. Этот файл описывает покрытие требований, а не заменяет критерии приёмки. Основание первоначального аудита: исходники `phz`, примеры, тесты и отчёт версии 0.4.0. Предыдущая рабочая итерация дала прогресс: разреженный FEM, 157 тестов, 15 контрольных расчётов и воспроизведение большого 3D-примера. Формальное доказательство этими тестами не установлено.

| Пункт | Текущее покрытие и необходимое доказательство приёмки |
|---|---|
| 0.1–0.3 | Есть локальный цикл постановка–дискретизация–расчёт–численные проверки–архив. Универсальность, публикационная и промышленная пригодность, формальная корректность всех результатов не доказаны. |
| 0.4 | Открытая задача ещё не выбрана и не решена. Нужны конкретный научный результат, поиск предшествующих работ и формальная проверка. |
| 1.1 | Есть Lark-парсер, типизированный ASCII DSL/IR для скалярных, векторных и тензорных систем и глобальный интеграл `intg`. В 0.7.0 векторный подкласс несжимаемого потока соединён с численной P2/P1-сборкой из слабой формы. Общие интегральные ядра, дробные производные и вычисление произвольных тензорных систем ещё не выполнены. |
| 1.2 | Автоматически выводятся объёмные и граничные члены Навье–Стокса в 2D/3D, включая симметричные напряжения. В 0.7.0 полученная форма используется для сборки смешанного остатка и якобиана. Проверены компонентные формулы, тождества Грина и численные изготовленные решения. Lean/Coq-сертификата нет. |
| 1.3 | Запрошенные пользователем массы проверяются численно. Автоматического поиска симметрий Ли/Нётер и генератора бесконечной иерархии КдФ нет. |
| 1.4 | В 0.5.0 добавлены точные рациональные таблицы RK, нелинейный порядок по корневым деревьям, достаточные полиномиальные свидетельства A-/L-устойчивости и вывод ограничения шага для косинусной теплопроводности. Для FEM условие пока зависит от несертифицированного обобщённого спектра. Есть численные исследования сходимости. Общие апостериорные оценки, Lean/Coq и полноценная автоматическая теорема сходимости не выполнены. |
| 2.1 | Косинусный базис, скалярный P1 и смешанный P2/P1. Библиотечный класс CPn строит непрерывные симплексные P1–P4, но полное пользовательское переключение всех решателей на P2–P4 ещё не выполнено. Нет остальных спектральных базисов, Hermite, Nedelec, Raviart–Thomas, вейвлетов и спектральных элементов. |
| 2.2 | Автовыбор базиса, DG и обученная модель отсутствуют. Нужны проверки выбора на гладкой и разрывной задачах. |
| 2.3 | Гаусс и отдельное повышение квадратуры доступны. В 0.8.0 добавлен автоматический подбор порядка по изменениям решения и действия операторов, с журналом и явным отказом подтверждения на пределе. Оценка эмпирическая, без гарантии для пропущенных правилом особенностей. Лобатто, MC/QMC, локальное адаптивное разбиение и интервальная оценка остаются невыполненными. |
| 3.1 | SciPy BDF/Radau/RK45/DOP853, собственный BE и фиксированные FE/RK2/RK4/SSPRK3. В 0.7.0 BE расширен на связанную смешанную DAE скорости и давления с разреженным Ньютоном. Временное уточнение поддерживает фиксированные методы. Нет IMEX/ETD/сохраняющих схем и полной адаптивности. Критерий жёсткой задачи без ручной настройки требует отдельного набора проверок. |
| 3.2 | Есть равномерное вложенное сгущение. Нет локального AMR, апостериорного управления и доказанного выигрыша в 10 раз на пограничном слое. |
| 3.3 | Вычисления локальные. MPI, GPU, декомпозиция и масштабирование до 1000 ядер не реализованы и не проверены. |
| 4.1 | Плотный и разреженный демпфированный Ньютон, заданные массы, локальный спектр. Нет Newton–Krylov, псевдодугового продолжения и deflation. Нужны ветви нескольких устойчивых и неустойчивых равновесий. |
| 4.2 | Полного бифуркационного анализа нет. Лоренц, Хопф, складки, гомоклинические/гетероклинические орбиты и коразмерность 2+ остаются открытыми работами. |
| 4.3 | Байесовские обратные задачи, MCMC/HMC/VI, GP/PCE и Sobol отсутствуют. Нужны синтетические эксперименты с интервалами и проверкой покрытия. |
| 5.1 | FNO/DeepONet/PINN/transformers отсутствуют. Ошибка менее 1% на независимом тесте не подтверждена. |
| 5.2 | Гибридное управление ML и классическим расчётом отсутствует. Выигрыш 10+ раз при одинаковой точности не подтверждён. |
| 5.3 | SINDy/символьная регрессия отсутствуют. Восстановление Навье–Стокса и вязкости с ошибкой до 5% не проверено. |
| 6.1 | Lean/Coq-доказательства не генерируются; Lean toolchain в PATH не найден. Символьная тождественность и численные тесты не засчитываются как формальное доказательство. Нужны компиляция без `sorry` и проверка аксиом. |
| 6.2 | Гарантированных интервальных оценок ошибки с направленным округлением нет. |
| 7.1 | Есть legacy Tkinter, CLI и Python-модули. Публичный стабильный SDK, Jupyter-интеграция, React/backend и REST/gRPC ещё не выполнены. |
| 7.2 | SLURM/PBS/Kubernetes и AWS/GCP/Azure из браузера отсутствуют. Проверка реального запуска требует соответствующей среды. |
| 7.3 | Версионируемый API плагинов и разрешение зависимостей отсутствуют. |
| 7.4 | Есть снимок Python-кода, среды, данных, сетки и контрольные суммы; локальный replay проверен. Нет Docker/Singularity и RO-Crate/PROV. Воспроизводимость через 10/20 лет сейчас не может быть подтверждена наблюдением. |
| 8.1 | Есть локальные примеры и верификация. Публичный сравнительный набор вплоть до турбулентности не опубликован. |
| 8.2 | Статьи не подготовлены и не опубликованы; приёмка требует реальных публикаций 3–5 работ, а не текстовых заготовок. |
| 8.3 | Открытая научная задача и формально проверенный новый результат отсутствуют. Сохраняется полный объём требования. |

## Последовательность развития

1. Проверяемое символьное ядро: анализ схем, математический DSL/IR, слабые формы и сертификаты. Затем подключение Lean/Coq и интервальных оценок с явно заданными предпосылками.
2. Базы и квадратуры, сохраняющие/IMEX/ETD-интеграторы, локальная адаптация и масштабируемая линейная алгебра.
3. Продолжение ветвей, deflation, бифуркации, обратные задачи и UQ.
4. ML-операторы, гибридная проверка и восстановление уравнений на воспроизводимых наборах данных.
5. Интерфейсы, плагины, контейнеры и provenance, кластерные среды и масштабирование.
6. Публичные сравнительные эксперименты, конкретная открытая задача, формальные результаты и публикации.

Пункты могут развиваться параллельно по архитектурным зависимостям, но ни один не исключается из цели. Завершение локального этапа не означает завершение полного ТЗ.

## Проверка этапа 0.5.0

Общий прогон: 198 тестов прошли с предупреждениями, превращёнными в ошибки. Отдельно проверены порядок на нелинейной ODE, отличие линейного и нелинейного порядка, чтение длинных десятичных коэффициентов без округления, полюса внутренних стадий, точное ограничение шага в 1D/2D/3D, нарушение ограничения явным Эйлером, зависимые от времени границы, временная сходимость, CLI, архив и replay. Исходники ядра и тестов проверены на ASCII и отсутствие комментариев; новые вычислительные функции имеют короткие имена с венгерскими префиксами.

Контрольный прогон: 16 расчётов, 5 исследований сходимости, 7 символьных профилей схем. На `rkheat` наблюдаемый порядок RK4 равен 4.06135. Для FE автоматически выведено `dt <= 8/(9*pi**2)`, а шаг 0.1 отмечен как нарушение. Все пять проверок повторения сохранённым кодом дали нулевую максимальную разность на этой машине. Подробные результаты: [VALIDATION.md](../generated/VALIDATION.md), машинные данные: [validation.json](../generated/validation.json).

Точные коэффициенты и полиномиальные тождества — подготовленный вход для дальнейшей формализации, но не Lean-доказательство. Полная цель остаётся незавершённой. Следующие работы по фундаменту: проверяемое представление PDE/слабых форм, перенос утверждений в proof assistant, апостериорные и интервальные оценки. Пункты исследования ветвей, AMR, ML, HPC, публикаций и открытой задачи сохраняются в полном объёме.

## Проверка этапа 0.6.0

Добавлены типизированный ASCII DSL/IR, автоматические слабые формы и граничные члены, скалярный численный мост и глобальные интегральные члены. Описание и точные границы поддержки: [DSL.md](DSL.md). Навье–Стокс выводится в 2D/3D как через векторный лапласиан, так и через симметричный тензор напряжений. Его численного решателя пока нет.

Общий прогон: 259 тестов прошли с предупреждениями, превращёнными в ошибки; провалов, ошибок и пропусков нет. Включены 61 новая проверка: компонентные формулы Навье–Стокса, точные тождества Грина на полиномах, переменные коэффициенты, приоритет операций, типы и недопустимый синтаксис, нелинейные интегралы, нелокальный якобиан, сохранение массы, аналитическое затухание моды, CLI, архив и replay. Зависимости согласованы; ASCII и отсутствие комментариев/docstring проверены для ядра и тестов.

Контрольный прогон: 18 расчётов, 6 исследований сходимости, 7 символьных профилей схем. Максимальная L2-ошибка изготовленного решения нелинейной диффузии `dsl2` равна `6.417667343331008e-10`. Дрейф массы объёмной нелокальной модели `int3` равен `1.3322676295501878e-15`. При сгущении её сетки с 64 до 343 узлов конечная L2-ошибка уменьшилась с `0.01680439999784665` до `0.007606373953579351`; квадратура оставалась фиксированной. Это наблюдение чувствительности к сетке, не гарантированная оценка ошибки.

Все шесть повторений сохранёнными копиями кода дали нулевую максимальную разность на этой машине. Новые файлы `weak.json` и `weak.md` входят в проверяемый архив. Полный отчёт: [VALIDATION.md](../generated/VALIDATION.md); машинные данные: [validation.json](../generated/validation.json).

Следующий фундаментальный этап: переход от скалярного моста к сборке смешанных пространств для векторных систем, оценка ошибки и формальные сертификаты при явно заданных предпосылках. Остальная матрица полного ТЗ остаётся действующей; цель не завершена.

## Проверка этапа 0.7.0

Добавлен численный путь Стокса и Навье–Стокса в 2D/3D: непрерывные пространства P2/P1, разреженная сборка остатка и его якобиана непосредственно из типизированной слабой формы, стационарный Ньютон и связанный Backward Euler для скорости и давления. Нелинейные значения и градиенты вычисляются в квадратурных точках каждого элемента. Поддерживаются полные границы Дирихле для скорости и нулевое среднее давление. Подробная постановка, команды и ограничения: [FLOW.md](FLOW.md).

Общий набор дал 299 пройденных тестов с предупреждениями как ошибками. После последних изменений повторены все 37 тестов смешанного пути; ошибок, провалов и пропусков нет. Проверены P1–P4, точные элементные интегралы, ориентация общих рёбер, изготовленные нелинейные решения в 2D/3D с нулевым начальным приближением, якобиан независимым численным дифференцированием, давление, меняющиеся во времени границы, несовместимый поток, начальная проекция скорости, энергия, CLI и архив. Сохранение архива при временной блокировке файлов в Windows отдельно проверено. Ядро и тесты используют ASCII без комментариев и docstring; зависимости согласованы.

Контрольный прогон версии 0.7.0 завершён: 22 расчёта, 9 исследований уточнения и 7 символьных профилей схем. Все восемь повторений сохранёнными копиями кода дали нулевую максимальную разность на этой машине, включая стационарный и нестационарный 3D-поток. Настольное окно проверено загрузкой смешанной модели, выводом слабой формы, решением и переключением на давление.

На неполиномиальном 3D-примере использованы 402, 2312 и 15468 коэффициентов. L2-ошибка скорости уменьшилась `0.493371 → 0.143181 → 0.0403313`, давления — `1.96582 → 0.663447 → 0.169180`, L2-дивергенция — `2.22565 → 1.39167 → 0.763443`. Квадратура и допуск Ньютона фиксированы. Слабые моменты несжимаемости имеют другую роль и сохраняются отдельно: они не заменяют точечную дивергенцию и ошибку поля. Наблюдаемые временные порядки BE по разностям решений равны `0.863189` для скорости и `1.00095` для давления; это численное наблюдение на данной серии, не теорема сходимости.

Архив дополнен пространствами скорости и давления, картами степеней свободы, матрицей дивергенции, последним слабым якобианом и историей Ньютона. Начальное состояние нестационарной задачи явно отмечено как не имеющее невязки шага BE. Повышение квадратуры, сгущение сетки и уменьшение шага времени исследуются раздельно; повышение порядка квадратуры не получает ложного порядка по формуле для деления шага пополам.

Полные числа и пути архивов: [VALIDATION.md](../generated/VALIDATION.md), [validation.json](../generated/validation.json). Открытые границы, локальное AMR, прочие базисы из ТЗ, MPI/GPU, анализ устойчивости смешанной DAE и Lean/Coq-сертификаты ещё не реализованы. Цель остаётся активной; завершён только описанный вычислительный этап.

## Проверка этапа 0.8.0

Добавлен автоматический подбор порядка квадратуры для скалярного FEM, косинусного пространства и смешанного потока. Он сравнивает поля и действие операторов на одних состояниях, проверяет начальные приближения, энергию и алгебраические поля, когда они присутствуют. Требуются последовательные подтверждения; исчерпание порядка без них сохраняет результат с `fail`, а не выдаёт подтверждение точности. Сетка и настройки решателя в ходе подбора фиксированы. Подробности: [QUADRATURE.md](QUADRATURE.md).

Общий прогон дал 333 пройденных теста без ошибок, провалов и пропусков. После расширения допустимого порядка до 20 отдельно прошли 35 целевых тестов, включая дополнительный тест высокого порядка; эти числа не складываются, поскольку наборы пересекаются. Ядро и тесты проверены на ASCII и отсутствие комментариев/docstring, установлена версия 0.8.0, зависимости согласованы.

Повторный общий контроль версии 0.8.0 также завершён: 22 расчёта, 9 исследований уточнения и 7 профилей схем. Вместе с шестью новыми расчётами автоподбора выполнены 28 контрольных расчётов. Итоговые поля четырёх повторённых архивов автоподбора совпали с нулевой максимальной разностью на этой машине; межплатформенная побитовая идентичность из этого не следует.

Устранены два найденных численных сбоя. Для P1 на тетраэдрах требуется квадратура как минимум третьего порядка, чтобы квадратурная матрица массы совпадала с используемой согласованной матрицей. Для BDF/Radau с алгебраическим Пуассоном используется вычисление якобиана с подбором допустимых возмущений; исходная проверка нейтральности не ослаблена, не нейтральные данные по-прежнему отклоняются.

Пройдены шесть контрольных расчётов автоподбора. Для экспоненциальной нагрузки FEM выбран порядок 10, для косинусного примера — 12, для нелинейного потока в 2D — 9, в 3D — 11. Для 3D сохранены прежние допуски и сетка; расширен предел квадратуры. Архив хранит фактическую модель, исходный запрос подбора и историю сравнений; повторение запускает весь подбор заново. Воспроизведение проверено для 3D-FEM, 3D-потока, косинусного примера и нестационарного FEM с численным допуском `1e-10`. Промежуточные контрольные суммы могут отличаться из-за округления между процессами и не подменяют сравнение чисел.

Настольное окно проверено для равновесия с автоподбором, расчёта во времени и переключения смешанного поля на давление. Краткая инструкция для локального запуска: [START.md](../START.md). Отчёт новых контролей: [AQ_VALIDATION.md](../generated/AQ_VALIDATION.md), машинные данные: [aq_validation.json](../generated/aq_validation.json).

Это контроль чувствительности к квадратуре, не гарантированный сертификат ошибки или замена сгущению сетки. Полное расширенное ТЗ остаётся незавершённым. По просьбе пользователя дальнейшая работа над полной целью приостановлена; текущая версия доведена для локального запуска.
)PHZDOC"},
    };
    return d;
}

J getModel(const std::string& s) {
    if (examples().count(s))
        return J::parse(examples().at(s));
    J j = load(s);
    if (j.contains("dom") && j["dom"].contains("path"))
        j["dom"]["path"] =
            (fs::absolute(fs::path(s)).parent_path() / j["dom"]["path"].get<std::string>())
                .string();
    return j;
}
}
