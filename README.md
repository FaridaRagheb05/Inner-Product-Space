# Inner-Product-Space
# Inner Product Verification Tool (C++)

## Overview

This project is a C++ program that **empirically tests whether a given function behaves like an inner product on \( \mathbb{R}^n \)**.

It checks the three defining axioms of an inner product:

1. **Positive Definiteness**  
   \( \langle v, v \rangle \ge 0 \) for all vectors \( v \), and equals 0 iff \( v = 0 \)

2. **Symmetry**  
   \( \langle u, v \rangle = \langle v, u \rangle \)

3. **Linearity (in first argument)**  
   \( \langle au + bv, w \rangle = a\langle u, w \rangle + b\langle v, w \rangle \)

> Important: This program does NOT prove a function is an inner product. It only performs randomized numerical testing. Passing all tests is strong evidence, but not a formal proof.

---

## Features

- Tests multiple candidate inner product functions
- Randomized vector generation for robustness
- Supports multiple dimensions (\( \mathbb{R}^n \))
- Built-in examples:
  - Standard dot product (valid)
  - Weighted dot product (valid if weights > 0)
  - Antisymmetric form (invalid)
  - Negative-weight product (invalid)
  - Nonlinear function (invalid)

- Detailed console output for each axiom
- Final summary verdict

---

## How It Works

The program:

1. Generates random vectors in \( \mathbb{R}^n \)
2. Applies the chosen inner product function
3. Checks:
   - positivity of \( \langle v,v \rangle \)
   - symmetry between \( \langle u,v \rangle \) and \( \langle v,u \rangle \)
   - linearity using random scalars and vectors
4. Tracks pass/fail statistics for each axiom
5. Prints a final verdict

---

## Compilation

Use any C++17-compatible compiler:

```bash
g++ -std=c++17 main.cpp -o inner_product_test
````

---

## Usage

```bash
./inner_product_test [options]
```

### Options

| Option                    | Description                              |
| ------------------------- | ---------------------------------------- |
| `-d, --dimension DIM`     | Dimension of vector space (default: 3)   |
| `-f, --function TYPE`     | Function type to test                    |
| `-w, --weights W1 W2 ...` | Weights for weighted dot product         |
| `-t, --tolerance EPS`     | Floating-point tolerance (default: 1e-9) |
| `-h, --help`              | Show help message                        |

---

## Function Types

### 1. Standard Dot Product

```text
<u,v> = u₁v₁ + u₂v₂ + ... + uₙvₙ
```

Expected: VALID inner product

---

### 2. Weighted Dot Product

```text
<u,v> = Σ wᵢ uᵢ vᵢ
```

Valid only if all weights wᵢ > 0

---

### 3. Negative Weights (Invalid)

```text
<u,v> = -u₁v₁ + u₂v₂ + ...
```

Fails positive definiteness

---

### 4. Antisymmetric Function (Invalid)

```text
<u,v> = u₁v₂ - u₂v₁
```

Fails symmetry

---

### 5. Nonlinear Function (Invalid)

```text
<u,v> = (u₁v₁)² + u₂v₂ + ...
```

Fails linearity

---

## Example Runs

### Test standard dot product

```bash
./inner_product_test -d 3 -f dot
```

### Test weighted inner product

```bash
./inner_product_test -d 4 -f weighted -w 1 2 3 4
```

### Test nonlinear function

```bash
./inner_product_test -d 2 -f nonlinear
```

---

## Output

For each axiom, the program prints:

* Individual test cases
* Computed values
* Pass/fail result per test
* Axiom-level summary

Final output:

```
VERDICT: VALID INNER PRODUCT
```

or

```
VERDICT: NOT A VALID INNER PRODUCT
```

---

## Limitations

* Randomized testing -> not a formal proof
* Finite sampling only
* Floating-point precision may affect borderline cases
* Cannot guarantee correctness for all vectors in ℝⁿ

---

## Educational Purpose

This tool is designed to help students:

* Understand inner product axioms
* See why certain functions fail algebraic properties
* Connect theory (linear algebra) with computational testing
* Debug intuition about vector space structure

---

## Author Notes

This project is intentionally built as a **verification simulator**, not a symbolic proof system. Its goal is to make abstract linear algebra concepts more concrete through computation.

```
```
