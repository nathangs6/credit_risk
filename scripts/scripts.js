// Note: must have mathjs loaded

/////////////////////
/// MATH ROUTINES ///
/////////////////////
function n(x) {
    return Math.exp(-(x**2) / 2) / Math.sqrt(2 * Math.PI)
}

function N(x) {
    return 0.5 * (1 - math.erf((-x)/(Math.sqrt(2) * 1)))
}

function central_difference(f, x0, dt) {
    return (f(x0+dt) - f(x0-dt))/2/dt
}

//////////////
/// LINALG ///
//////////////
class Vec2 {
    constructor(x0, x1) {
        this.x0 = x0
        this.x1 = x1
    }

    scalarAdd(a) {
        return Vec2(this.x0 + a, this.x1 + a)
    }

    vecAdd(v) {
        return Vec2(this.x0 + v.x0, this.x1 + v.x1)
    }

    scalarMul(a) {
        return Vec2(a * this.x0, a * this.x1)
    }

    norm() {
        return Math.sqrt(this.x0**2 + this.x1**2)
    }
}

class Mat22 {
    constructor(x00, x01, x10, x11) {
        this.x00 = x00
        this.x01 = x01
        this.x10 = x10
        this.x11 = x11
    }

    det() {
        return this.x00 * this.x11 - this.x01 * this.x10
    }

    scalarMul(a) {
        return Mat22(a*this.x00, a*this.x01, a*this.x10, a*this.x11)
    }

    vecMul(v) {
        return Vec2(this.x00 * v.x0 + this.x01 * v.x1, this.x10 * v.x0 + this.x11 * v.x1)
    }

    matMul(M) {
        m00 = this.x00 * M.x00 + this.x01 * M.x10
        m01 = this.x00 * M.x01 + this.x01 * M.x11
        m10 = this.x10 * M.x00 + this.x11 * M.x10
        m11 = this.x10 * M.x01 + this.x11 * M.x11
        return Mat22(m00, m01, m10, m11) 
    }

    inverse() {
        d = this.det()
        return Mat22(this.x11, this.x01, this.x10, this.x11).scalarMul(-d)
    }

    transpose() {
        return Mat22(this.x00, this.x10, this.x01, this.x11)
    }
}


////////////////////////
/// OPTION FUNCTIONS ///
////////////////////////
function getStrike(K, r, t) {
    return K * Math.exp(-r*t)
}
function dVals(V, A, v, t) {
    const a = Math.log(V/A)
    const b = 0.5 * v**2 * t
    const c = v * Math.sqrt(t)
    const d1 = (a + b) / c
    return [d1, d1 - v * Math.sqrt(t)]
}

function priceOption(V, A, v, t) {
    const [d1, d2] = dVals(V, A, v, t)
    return V * N(d1) - A * N(d2)
}

function delta(V, A, v, t) {
    const [d1, d2] = dVals(V, A, v, t)
    return N(d1)
}

/////////////////////////////
/// CREDIT RISK FUNCTIONS ///
/////////////////////////////
function compare_sigmaV(V, vV, S, vS, delt) {
    return vS * S - vV * V * delt
}

function compare_S(S, S0) {
    return S - S0
}

function computeJacobian(V, vV, A, t, vS) {
    const [d1, d2] = dVals(V, A, V, t)
    const nd1 = n(d1)
    const Nd1 = N(d1)

    SdV = Nd1 // delta
    Sdv = V * nd1 * Math.sqrt(t) // vega

    gdV = (vS - vV) * Nd1 - nd1/Math.sqrt(t)
    gdv = vS * Sdv - V * (Nd1 + nd1*d2)

    return math.matrix([[SdV, Sdv], [gdV, gdv]])
}

function findX(s, v, w, D) {
    let sv = math.multiply(v, s)
    let wmsv = math.add(w, math.multiply(sv, -1))

    let a = math.norm(wmsv)**2
    let b = 2 * math.dot(sv, wmsv)
    let c = math.norm(sv)**2 - D**2

    return (-b + Math.sqrt(b**2 - 4*a*c))/(2*a)
}

function find_pdl_step(D, V, vV, A, t, vS, S0) {
    let J = computeJacobian(V, vV, A, t, vS)
    let JT = math.transpose(J)
    let S = priceOption(V, A, vV, t)
    let delt = delta(V, A, vV, t)
    let f1 = compare_S(S, S0)
    let f2 = compare_sigmaV(V, vV, S, vS, delt)
    let nfx = math.matrix([-f1, -f2])

    let d_gn = math.multiply(math.multiply(math.inv(math.multiply(JT, J)), JT), nfx)
    let n_gn = math.norm(d_gn)
    let d_sd = math.multiply(JT, nfx)
    let n_sd = math.norm(d_sd)
    let dk

    if (n_gn < D) {
        dk = d_gn
    } else {
        let s = n_sd**2 / math.norm(math.multiply(J, d_sd))**2
        let D2 = s * n_sd
        if (n_gn > D2 && n_sd > D2) {
            dk = math.multiply(d_sd, D / n_sd)
        } else {
            dk = findX(d_sd, d_gn, s, D)
        }
    }
    return dk.toArray()
}

function solve(V0, v0, A, t, S0, vS, max_iter, D) {
    let V = V0
    let vV = v0
    let dV
    let dvV
    for (let i = 0; i < max_iter; i++) {
        [dV, dvV] = find_pdl_step(D, V, vV, A, t, vS, S0)
        V = V + dV
        vV = vV + dvV
    }
    return [V, vV]
}

function findSolutions() {
    let S0 = document.getElementById("equity").value
    let vS = document.getElementById("equity-volatility").value
    let K = document.getElementById("debt").value
    let r = document.getElementById("risk-free-rate").value
    S0 = Number(S0)
    vS = Number(vS) / 100
    K = Number(K)
    r = Number(r) / 100
    let t = 1
    let A = getStrike(K, r, t)
    let V0 = S0 + K
    let v0 = vS

    let max_iter = 1000
    let D = 1

    let [V, vV] = solve(V0, v0, A, t, S0, vS, max_iter, D)
    document.getElementById("assets").textContent = math.round(V, 4)
    document.getElementById("asset-volatility").textContent = math.round(vV * 100, 4)
    document.getElementById("equity-output").textContent = math.round(priceOption(V, A, vV, t), 4)
    document.getElementById("delta").textContent = math.round(delta(V, A, vV, t), 4)
}