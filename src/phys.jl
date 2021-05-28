# Consider using Unitful.jl.

struct PhysUnit
    # Quantities whose units are fundamental constants
    ε::Float  # permittivity unit
    μ::Float  # permeability unit
    c::Float  # speed unit
    η::Float  # impedance unit

    # General quantities (not specific to electromagnetism)
    L::Float  # length unit
    ω::Float  # angular frequency unit

    # Electromagnetic quantities, in order of frequency of appearance in Maxwell's equations
    E::Float  # E-field unit
    H::Float  # H-field unit

    J::Float  # electric current density unit
    M::Float  # magnetic current density unit

    S::Float  # Poynting vector unit
    P::Float  # power unit

    V::Float  # voltage unit
    I::Float  # electric current unit

    function PhysUnit(L₀::Real)
        # The length unit L₀ is set by users.
        L₀ > 0 || throw(ArgumentError("L₀ = $L₀ must be positive."));

        # The E-field unit is arbitrarily chosen as E₀ = 1 V/m.
        E₀ = 1.0  # E-field is measured in V/m

        # The following quantities are measured in their fundamental constants:
        # - Permittivity is measured in vacuum permittivity ε₀
        # - Permeability is measured in vacuum permeability μ₀
        # - Speed is measured in speed c₀ of light in vacuum
        # - Impedance is measured in vacuum impedance η₀ = √(μ₀/ε₀)

        # The units of all the other quantities are derived from the above units, such that
        # they are consistent with the above units in Maxwell's equations.  Among them, the
        # following units are independent of the length unit L₀:
        H₀ = E₀ / η₀  # H-field is measured in in E₀/η₀ (L₀-independent)
        S₀ = E₀ * H₀  # Poynting vector is measured in E₀H₀ (L₀-independent)

        # The following units are dependent on the length unit L₀:
        ω₀ = c₀ / L₀  # frequency is measured in c₀/L₀
        J₀ = H₀ / L₀  # electric current density is measured in H₀/L₀
        M₀ = E₀ / L₀  # magnetic current density is measured in E₀/L₀
        I₀ = J₀ * L₀^2  # electric current is measured in J₀L₀²
        V₀ = E₀ * L₀  # voltage is measured in E₀L₀
        P₀ = I₀ * V₀  # power is measured in I₀V₀ = S₀L₀²

        new(ε₀, μ₀, c₀, η₀, L₀, ω₀, E₀, H₀, J₀, M₀, S₀, P₀, V₀, I₀)
    end
end
