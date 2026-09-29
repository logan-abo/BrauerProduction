#############################################################################
# MullerTreeTemplate.g
#
# Constructs the Brauer trees for faithful weight-1 blocks of
# the Schur cover of A_n, using the trees in Muller's Tables 4, 5, and 6.
#
# For epsilon = +1, plus and minus are lists of character positions and zero
# is the pair of exceptional character positions.
# For epsilon = -1, plus and minus are lists of split-character position pairs
# and zero is a single character position.
#############################################################################


#############################################################################
# Table 4
#
# epsilon = +1
#
#                  lambda_0,1/2
#              /                  \
#       lambda_t^- ... lambda_1^-   lambda_s^+ ... lambda_1^+
#
# The two characters corresponding to lambda_0 are exceptional.
#############################################################################

MullerTreePositive := function(plus, zero, minus, name)
    return rec(
        name := name,

        center := zero,

        branches := [
            Reversed(minus),
            Reversed(plus)
        ],

        exceptionalMultiplicity := 2
    );
end;


#############################################################################
# Tables 5 and 6
#
# epsilon = -1
#
# The center is lambda_0^(0).
# Each lambda_i^+ and lambda_i^- splits into two characters upon restriction
# to A~_n.
#
# The distinction between Tables 5 and 6 is only the geometric placement
# of the four branches; as an abstract labelled tree, the following four
# branches give the required tree.
#############################################################################

MullerTreeNegative := function(plus, zero, minus, name)
    local plus1, plus2, minus1, minus2;

    plus1 := List(
        plus,
        x -> x[1]
    );

    plus2 := List(
        plus,
        x -> x[2]
    );

    minus1 := List(
        minus,
        x -> x[1]
    );

    minus2 := List(
        minus,
        x -> x[2]
    );

    return rec(
        name := name,

        center := [zero],

        branches := [
            Reversed(plus1),
            Reversed(minus1),
            Reversed(plus2),
            Reversed(minus2)
        ],

        exceptionalMultiplicity := 1
    );
end;


#############################################################################
# Main Müller tree template.
#
# epsilon = +1:
#     Table 4
#
# epsilon = -1:
#     Tables 5/6
#
#############################################################################

MullerTreeTemplate := function(plus, zero, minus, epsilon, name)

    if epsilon = 1 then

        return MullerTreePositive(
            plus,
            zero,
            minus,
            name
        );

    elif epsilon = -1 then

        return MullerTreeNegative(
            plus,
            zero,
            minus,
            name
        );

    else

        Error(
            "MullerTreeTemplate: epsilon must be +1 or -1"
        );

    fi;

end;