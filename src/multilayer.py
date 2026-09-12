"""Arbitrary stacks of layers on a semi-infinite substrate.

`forward_model.Sample` describes one film on a substrate, which covers the
simplest case. Real samples rarely stop there: a metal transducer deposited
for pump-probe measurement, a native oxide, an adhesion layer, or a damaged
region near the surface all add layers, each with its own contact resistance.

A stack is described from the heated front face inwards:

    front surface (optional losses)
    layer 0
    resistance 0
    layer 1
    resistance 1
    ...
    layer n-1
    resistance n-1
    substrate (semi-infinite, entering through its impedance)

There must be exactly one resistance per layer, the last one sitting between
the deepest layer and the substrate. Set any of them to zero for perfect
contact.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

import forward_model as fm
import laplace as lp
import quadrupoles as q

__all__ = [
    "Layer",
    "Stack",
    "stack_matrix",
    "response",
    "modulated_response",
    "pulsed_response",
]


@dataclass
class Layer:
    """One layer of the stack.

    Setting `lam_back` and `rho_c_back` makes the layer graded, with a linear
    metaproperty profile between the front and back values.
    """

    lam: float
    rho_c: float
    thickness: float
    lam_back: float | None = None
    rho_c_back: float | None = None

    @property
    def a(self) -> float:
        return q.diffusivity(self.lam, self.rho_c)

    @property
    def b(self) -> float:
        return q.effusivity(self.lam, self.rho_c)

    @property
    def b_back(self) -> float | None:
        if self.lam_back is None or self.rho_c_back is None:
            return None
        return q.effusivity(self.lam_back, self.rho_c_back)

    @property
    def xi(self) -> float:
        return q.xi_from_thickness(self.thickness, self.a)

    @property
    def diffusion_time(self) -> float:
        return self.thickness ** 2 / self.a

    @property
    def is_graded(self) -> bool:
        return self.b_back is not None

    def matrix(self, p) -> np.ndarray:
        if self.is_graded:
            return q.graded_linear_layer(p, self.b, self.b_back, self.xi, form="T")
        return q.homogeneous_wall(p, self.lam, self.rho_c, self.thickness)


@dataclass
class Stack:
    """A sequence of layers closed by a semi-infinite substrate."""

    layers: list[Layer]
    sub_lam: float = 35.0
    sub_rho_c: float = 3.03e6
    contact_resistances: list[float] = field(default_factory=list)
    front_losses: float = 0.0

    def __post_init__(self):
        if not self.layers:
            raise ValueError("a stack needs at least one layer")
        if not self.contact_resistances:
            self.contact_resistances = [0.0] * len(self.layers)
        if len(self.contact_resistances) != len(self.layers):
            raise ValueError(
                "there must be exactly one contact resistance per layer, "
                "the last one lying between the deepest layer and the substrate"
            )

    @property
    def sub_b(self) -> float:
        return q.effusivity(self.sub_lam, self.sub_rho_c)

    @property
    def total_thickness(self) -> float:
        return sum(l.thickness for l in self.layers)

    @property
    def total_diffusion_time(self) -> float:
        """Square of the total Liouville thickness.

        The Liouville thicknesses add, not the diffusion times, which is why
        this is the square of a sum and not a sum of squares.
        """
        return sum(l.xi for l in self.layers) ** 2

    @property
    def characteristic_frequency(self) -> float:
        return 1.0 / (2.0 * np.pi * self.total_diffusion_time)


# --------------------------------------------------------------------------
# Forward model
# --------------------------------------------------------------------------

def stack_matrix(p, s: Stack) -> np.ndarray:
    """Quadrupole of the whole stack, front face to substrate."""
    m = None
    for layer, r in zip(s.layers, s.contact_resistances):
        block = layer.matrix(p)
        if r != 0.0:
            block = block @ fm.interface_resistance(p, r)
        m = block if m is None else m @ block
    return m


def response(p, s: Stack, power=1.0):
    """Front-face temperature in the Laplace domain, per unit area."""
    m = stack_matrix(p, s)
    z = q.semi_infinite_impedance(p, s.sub_b)
    return q.front_face_temperature(m, z, power=power, h=s.front_losses)


def modulated_response(freq, s: Stack, power=1.0):
    freq = np.asarray(freq, dtype=float)
    theta = response(2j * np.pi * freq, s, power=power)
    return np.abs(theta), np.degrees(np.angle(theta))


def pulsed_response(t, s: Stack, energy=1.0, n_stehfest: int = 12):
    return lp.stehfest_inverse(
        lambda pp: np.real(response(pp, s, power=energy)), t, n=n_stehfest
    )
