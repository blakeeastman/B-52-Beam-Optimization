use duplicate::duplicate;
use strum_macros::EnumIter;

//Loads
const VIN: f64 = 800.0;
const VOUT: f64 = 1400.0;
const MDOT: f64 = 30.0;
const ENG_MASS: f64 = 16000.0;
const ENG1_LOC: f64 = 501.0;
const ENG2_LOC: f64 = 879.0;

const WEIGHT_FUEL: f64 = 252.0;
const FORCE_LIFT: f64 = 720.0;

// materials
// [density (lb/in^3), yield strength (psi), elastic modulus (psi), cost ($/lb), SigFb (ksi), A (#), B (#)]
const STEEL_1018: [f64; 7] = [
    0.284,
    76000.0,
    30000000.0,
    1.09,
    105.0,
    97.9684641113648,
    -0.1,
];
const STAINLESS_STEEL_17_4PH: [f64; 7] = [
    0.286,
    165000.0,
    28000000.0,
    4.65,
    257.474868044485,
    257.0,
    -0.095,
];
const STEEL_SAE_4340: [f64; 7] = [0.283, 132000.0, 29000000.0, 1.22, 237.0, 238.0, -0.0977];
const ALUMINUM_7075_T6: [f64; 7] = [
    0.102,
    73000.0,
    10400000.0,
    9.56,
    108.0,
    192.900038382796,
    -0.143,
];
const ALUMINUM_2024_T4: [f64; 7] = [0.100, 47000.0, 10600000.0, 5.35, 91.5, 122.0, -0.102];
const ALUMINUM_6061_T6: [f64; 7] = [
    0.0975,
    40000.0,
    10000000.0,
    2.92,
    85.0,
    101.666666666667,
    -0.107,
];
const TITANIUM_ALLOY_TI_6AL_4V: [f64; 7] = [0.16, 128000.0, 16500000.0, 61.50, 249.0, 274.0, -1.04];

#[derive(Copy, Clone, EnumIter)]
pub enum Material {
    Steel1018,
    StainlessSteel174PH,
    SteelSae4340,
    Aluminum7075T6,
    Aluminum2024T4,
    Aluminum6061T6,
    TitaniumAlloyTi6AL4V,
}

fn get_material(material: Material) -> [f64; 7] {
    return match material {
        Material::Steel1018 => STEEL_1018,
        Material::StainlessSteel174PH => STAINLESS_STEEL_17_4PH,
        Material::SteelSae4340 => STEEL_SAE_4340,
        Material::Aluminum7075T6 => ALUMINUM_7075_T6,
        Material::Aluminum2024T4 => ALUMINUM_2024_T4,
        Material::Aluminum6061T6 => ALUMINUM_6061_T6,
        Material::TitaniumAlloyTi6AL4V => TITANIUM_ALLOY_TI_6AL_4V,
    };
}

//Beam Phenotype Structs
#[derive(Copy, Clone)]
pub struct RectBeam {
    pub Material: Material,
    pub Length: f64,
    pub Width: f64,
    pub Height: f64,
    pub Thickness: f64,
}

#[derive(Copy, Clone)]
pub struct TBeam {
    pub Material: Material,
    pub Length: f64,
    pub Width: f64,
    pub Height: f64,
    pub StemThickness: f64,
    pub FlangeThickness: f64,
}

#[derive(Copy, Clone)]
pub struct IBeam {
    pub Material: Material,
    pub Length: f64,
    pub Width: f64,
    pub Height: f64,
    pub CenterThickness: f64,
    pub FlangeThickness: f64,
}

//Getting Beam Area
trait Area {
    fn area(&self) -> f64;
}

impl Area for RectBeam {
    fn area(&self) -> f64 {
        let inner_rect_area =
            (self.Width - 2.0 * self.Thickness) * (self.Height - 2.0 * self.Thickness);
        return (self.Width * self.Height) - inner_rect_area;
    }
}

impl Area for TBeam {
    fn area(&self) -> f64 {
        return self.Width * self.FlangeThickness
            + (self.Height - self.FlangeThickness) * self.StemThickness;
    }
}

impl Area for IBeam {
    fn area(&self) -> f64 {
        return 2.0 * self.Width * self.FlangeThickness
            + self.CenterThickness * (self.Height - 2.0 * self.FlangeThickness);
    }
}

//Getting Beam Weight
pub trait Weight {
    fn weight(&self) -> f64;
}

#[duplicate(beam_type; [RectBeam] ; [TBeam] ; [IBeam])]
impl Weight for beam_type {
    fn weight(&self) -> f64 {
        let volume = self.area() * self.Length;
        return volume * get_material(self.Material)[0];
    }
}

// Getting the cost of beams.
pub trait Cost {
    fn cost(&self) -> f64;
}

#[duplicate(beam_type; [RectBeam] ; [TBeam] ; [IBeam])]
impl Cost for beam_type {
    fn cost(&self) -> f64 {
        let weight = self.weight();
        return weight * get_material(self.Material)[3];
    }
}

// 2nd Moment of Area about x-axis.
trait Ix {
    fn ix(&self) -> f64;
}

impl Ix for RectBeam {
    fn ix(&self) -> f64 {
        let ix_big = self.Width * f64::powf(self.Height, 3.0) / 12.0;
        let ix_small = (self.Width - 2.0 * self.Thickness)
            * f64::powf(self.Height - 2.0 * self.Thickness, 3.0)
            / 12.0;
        return ix_big - ix_small;
    }
}

impl Ix for TBeam {
    fn ix(&self) -> f64 {
        let ix_stem =
            self.StemThickness * f64::powf(self.Height - self.FlangeThickness, 3.0) / 12.0;
        let ix_flange = self.Width * f64::powf(self.FlangeThickness, 3.0) / 12.0;
        let stem_displacement =
            self.y_bend() + (self.Height - self.FlangeThickness) / 2.0 + self.FlangeThickness;
        let flange_displacement = self.FlangeThickness / 2.0 + self.y_bend();
        let stem_displacement_ix = f64::powf(stem_displacement, 2.0)
            * (self.Height - self.FlangeThickness)
            * self.StemThickness;
        let flange_displacement_ix =
            f64::powf(flange_displacement, 2.0) * self.Width * self.FlangeThickness;
        return ix_stem + ix_flange + stem_displacement_ix + flange_displacement_ix;
    }
}

impl Ix for IBeam {
    fn ix(&self) -> f64 {
        let ix_small_rectangle = (self.Width - self.CenterThickness) / 2.0
            * f64::powf(self.Height - 2.0 * self.FlangeThickness, 3.0)
            / 12.0;
        let ix_large_rectangle = (self.Width) * f64::powf(self.Height, 3.0) / 12.0;
        return ix_large_rectangle - 2.0 * ix_small_rectangle;
    }
}

// Handle 2nd moment of area about y axis.
trait Iy {
    fn iy(&self) -> f64;
}

impl Iy for RectBeam {
    fn iy(&self) -> f64 {
        let iy_big = self.Height * f64::powf(self.Width, 3.0) / 12.0;
        let iy_small = (self.Height - 2.0 * self.Thickness)
            * f64::powf(self.Width - 2.0 * self.Thickness, 3.0)
            / 12.0;
        return iy_big - iy_small;
    }
}

impl Iy for TBeam {
    fn iy(&self) -> f64 {
        let iy_flange = f64::powf(self.Width, 3.0) * self.FlangeThickness / 12.0;
        let iy_stem =
            f64::powf(self.StemThickness, 3.0) * (self.Height - self.FlangeThickness) / 12.0;
        return iy_flange + iy_stem;
    }
}

impl Iy for IBeam {
    fn iy(&self) -> f64 {
        let iy_flange = self.FlangeThickness * f64::powf(self.Width, 3.0) / 12.0;
        let iy_stem = (self.Height - 2.0 * self.FlangeThickness)
            * f64::powf(self.CenterThickness, 3.0)
            / 12.0;
        return iy_stem + 2.0 * iy_flange;
    }
}

//Implement X Bend for beams

trait XBend {
    fn x_bend(&self) -> f64;
}

#[duplicate(beam_type; [RectBeam] ; [TBeam] ; [IBeam])]
impl XBend for beam_type {
    fn x_bend(&self) -> f64 {
        return self.Width / 2.0;
    }
}

//Implement Y bend for beams
trait YBend {
    fn y_bend(&self) -> f64;
}

#[duplicate(beam_type; [RectBeam]; [IBeam])]
impl YBend for beam_type {
    fn y_bend(&self) -> f64 {
        return -self.Height / 2.0;
    }
}

impl YBend for TBeam {
    fn y_bend(&self) -> f64 {
        let weighted_stem = (self.FlangeThickness + (self.Height - self.FlangeThickness) / 2.0)
            * self.StemThickness
            * (self.Height - self.FlangeThickness);
        let weighted_flange = (self.FlangeThickness / 2.0) * (self.FlangeThickness * self.Width);
        return -(weighted_stem + weighted_flange) / self.area();
    }
}

//Handle Stress:
pub trait Stress {
    fn stress_vert(&self) -> f64;
    fn stress_horz(&self) -> f64;
    fn total_stress(&self) -> f64;
    fn factor_of_safety(&self) -> f64;
    fn flight_hours(&self) -> f64;
    fn vertical_deflection(&self) -> f64;
}

#[duplicate(beam_type; [RectBeam] ; [TBeam] ; [IBeam])]
impl Stress for beam_type {
    fn stress_vert(&self) -> f64 {
        let m_eng1 = -ENG1_LOC * ENG_MASS;
        let m_eng2 = -ENG2_LOC * ENG_MASS;

        let m_lift = (FORCE_LIFT * self.Length) / 2.0 * self.Length / (3.0);
        let m_fuel = -(WEIGHT_FUEL * self.Length) / 2.0 * self.Length / 3.0;
        let m_weight = -(self.weight()) * self.Length / 2.0;

        return -(m_eng1 + m_eng2 + m_lift + m_fuel + m_weight) * self.y_bend() / self.ix();
    }
    fn stress_horz(&self) -> f64 {
        let thrust_eng = (VOUT - VIN) * MDOT;

        let m_eng1 = thrust_eng * ENG1_LOC;
        let m_eng2 = thrust_eng * ENG2_LOC;

        return (m_eng1 + m_eng2) * self.x_bend() / self.iy();
    }
    fn total_stress(&self) -> f64 {
        return self.stress_vert() + self.stress_horz();
    }
    fn factor_of_safety(&self) -> f64 {
        let yield_strength = get_material(self.Material)[1];
        return (yield_strength / self.total_stress());
    }
    fn flight_hours(&self) -> f64 {
        let material = get_material(self.Material);
        let g4 = material[4];
        let h4 = material[5];
        let i4 = material[6];
        let relevant_stress = self.stress_vert();
        let amplitude = relevant_stress / 2.0;
        let numerator = 2.33 * amplitude / 1000.0;
        let denominator = (1.0 - (2.33 * relevant_stress / 1000.0) / g4) * h4;
        return f64::powf(numerator / denominator, 1.0 / i4) / 27.0;
    }
    fn vertical_deflection(&self) -> f64 {
        let material = get_material(self.Material);
        let length_tofour = f64::powf(self.Length, 4.0);
        let ix = self.ix();
        let def_weight = -(self.weight() / self.Length) * length_tofour / 8.0 / material[2] / ix;
        let def_fuel = -WEIGHT_FUEL * length_tofour / 30.0 / material[2] / ix;
        let def_lift = FORCE_LIFT * length_tofour / 30.0 / material[2] / ix;
        let def_eng1 = -ENG_MASS * f64::powf(ENG1_LOC, 2.0) * (3.0 * self.Length - ENG1_LOC)
            / 6.0
            / material[2]
            / ix;
        let def_eng2 = -ENG_MASS * f64::powf(ENG2_LOC, 2.0) * (3.0 * self.Length - ENG2_LOC)
            / 6.0
            / material[2]
            / ix;
        return def_fuel + def_lift + def_eng1 + def_eng2 + def_weight;
    }
}
