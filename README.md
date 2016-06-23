# utvpi
A new way to generate integral solutions for Unit Two Variable Per Inequality constraint systems, in O(m*n) time.

Major theoretical work by Dr. K. Subramani and Piotr Wojciechowski. Implementation by Zachary Santer. Research conducted at West Virginia University, and funded by the National Science Foundation.

All C complies with the C11 specification and uses no platform-specific function calls (unless you define an __HPC__ macro when compiling subWojInt and lahiri). (Anonymous structs and unions are used in lahiri. Otherwise, C99 covers everything.)
