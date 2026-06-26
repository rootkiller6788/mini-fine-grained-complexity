# Gap Report: mini-equivalence-classes

## Priority 1 (High): Missing Core Knowledge

*None.* All L1-L8 knowledge areas are complete.

## Priority 2 (Medium): Partial Knowledge

| Item | Layer | Current State | Action Needed |
|------|-------|--------------|---------------|
| Fast Matrix Multiplication (omega < 2.372) | L5 | Fallback to naive min-plus | Implement Strassen-like min-plus with Z-strassen reduction |
| Gronlund-Pettie fully optimized 3SUM | L5 | Fallback to O(n^2) sort+two-pointer | Implement O(n^2/polylog n) with full hashing |

## Priority 3 (Low): Research Frontiers

| Item | Layer | Current State | Action Needed |
|------|-------|--------------|---------------|
| Fine-grained approximability | L9 | Documented only | Implement constant-factor approximation algorithms |
| Meta-complexity of fine-grained reductions | L9 | Not covered | Add theoretical analysis of reduction optimality |
| Quantum fine-grained complexity | L9 | Not covered | Add quantum algorithm lower bounds |

## Deferred Items
- Lean 4 formal proofs (separate module)
- Experimental validation against GPU/FPGA implementations