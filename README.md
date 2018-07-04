# Reciprocals

Trying to perform 64:32 bit division faster than native x64 DIV instruction.

Current approach:
- Get initial approximation of the reciprocal from LUT table
- If divisor is less than 2^31, do one Newton-Raphson iterations to get reciprocal with full 64-bit precision
- Multiply by reciprocal
- Correct the result

LUT table stores piecewise cubic polynomial approximation of 1/x. For each interval it stores 4 first coefficients of Taylor series for 1/x: C0+C1\*x+C2\*x^2+C3\*x^3. The origin point is taken in the middle of each inteval.
- C0 has 38 bits (37 bits stored explicitly)
- C1 has 27 bits
- C2 has 20 bits
- C3 has 12 bits