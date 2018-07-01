# Reciprocals

Trying to perform 64:32 bit division faster than native x64 DIV instruction.

Current approach:
- Get initial approximation of the reciprocal from LUT table
- Do one or two Newton-Raphson iterations to get reciprocal with full 64-bit precision
- Multiply by reciprocal
- Correct the result

LUT table stores piecewise linear approximation of 1/x. For each interval it stores:
- Value of 1/x in the middle point with 33 bits of precision (32 bits stored, MSB bit is always 1)
- Two linear coefficients used to approximate values between middle point and edges of the interval (17 bits of precision for each, 16 bits stored, MSB can be restored from interval index because these coefficients are strictly monotonic)
