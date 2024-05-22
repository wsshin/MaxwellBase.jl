# Consider using Unitful.jl.
export MaxwellUnit

struct MaxwellUnit
    # Quantities whose units are fundamental constants
    ε::Float64  # permittivity unit
    μ::Float64  # permeability unit
    c::Float64  # speed unit
    η::Float64  # impedance unit

    # General quantities (not specific to electromagnetism)
    L::Float64  # length unit
    ω::Float64  # angular frequency unit
    k::Float64  # wavenumber unit

    # Electromagnetic quantities, in order of frequency of appearance in Maxwell's equations
    E::Float64  # E-field unit
    H::Float64  # H-field unit

    J::Float64  # electric current density unit
    M::Float64  # magnetic current density unit

    S::Float64  # Poynting vector unit
    P::Float64  # power unit
    u::Float64  # energy density unit
    U::Float64  # energy unit

    V::Float64  # voltage unit
    I::Float64  # electric current unit

    function MaxwellUnit(L₀)
        # The length unit L₀ is set by users.
        @argcheck 1L₀ > 0  # L₀ must be positive

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
        u₀ = S₀ / c₀  # energy density is measured in ε₀E₀² = μ₀H₀² = E₀H₀/c₀ (L₀-independent)

        # The following units are dependent on the length unit L₀:
        ω₀ = c₀ / L₀  # frequency is measured in c₀/L₀
        k₀ = 1 / L₀  # wavenumber is measured in 1/L₀
        J₀ = H₀ / L₀  # electric current density is measured in H₀/L₀
        M₀ = E₀ / L₀  # magnetic current density is measured in E₀/L₀
        I₀ = J₀ * L₀^2  # electric current is measured in J₀L₀²
        V₀ = E₀ * L₀  # voltage is measured in E₀L₀
        P₀ = I₀ * V₀  # power is measured in I₀V₀ = S₀L₀²
        U₀ = u₀ * L₀^3  # energy is measured in u₀L₀³ = P₀L₀/c₀

        new(ε₀, μ₀, c₀, η₀, L₀, ω₀, k₀, E₀, H₀, J₀, M₀, S₀, P₀, u₀, U₀, V₀, I₀)
    end
end
