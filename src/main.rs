use rand::prelude::*;
use rsgenetic::pheno::*;
use rsgenetic::sim::select::*;
use rsgenetic::sim::seq::Simulator;
use rsgenetic::sim::*;
use std::string::String;
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread;
use strum::IntoEnumIterator;

mod beams;

//Modify the scoring function to change that the algorithm considers good. Recompile with a release version so it is optimized, run, and wait. It will spit out the best beam it finds.

use beams::*;
//GA specs
const POP_SIZE: u32 = 1000;
const POP_SURVIVORS: usize = 100;
const GENETIC_ITERS: u64 = 5000;

//dimensions
const LENGTH: f64 = 1257.0;
const HEIGHT: f64 = 34.0;

const WIDTH_MIN: f64 = 22.0;
const WIDTH_MAX: f64 = 38.0;

//requirements
const FOS_MAX: f64 = 2.2;
const FOS_MIN: f64 = 1.4;

const DEFLECTION_MAX: f64 = 70.0;

const PRICE_MAX: f64 = 500000.00;

const FATIGUE_FLIGHT_HOURS_MAX: f64 = 500000.0;
const FATIGUE_FLIGHT_HOURS_MIN: f64 = 42000.0;

const WEIGHT_MAX: f64 = 78000.0;

//Scoring design parameters.
fn score(
    length: f64,
    height: f64,
    width: f64,
    weight: f64,
    cost: f64,
    factor_of_safety: f64,
    deflection: f64,
    flight_hours: f64,
) -> i64 {
    let lengthcontrib = match length == LENGTH {
        true => 0 as i64,
        false => -i64::MAX / 10,
    };
    let heightcontrib = match height == HEIGHT {
        true => 0 as i64,
        false => -i64::MAX / 10,
    };
    let widthcontrib = match width < WIDTH_MAX && width > WIDTH_MIN {
        true => 0 as i64,
        false => -i64::MAX / 10,
    };
    let weightcontrib = match weight < WEIGHT_MAX {
        true => -weight as i64,
        false => -i64::MAX / 10,
    };
    let costcontrib = match cost < PRICE_MAX {
        true => (PRICE_MAX / 3.0 - cost / 3.0) as i64,
        false => -i64::MAX / 10,
    };
    let factorcontrib = match factor_of_safety < FOS_MAX && factor_of_safety > FOS_MIN {
        true => factor_of_safety*WEIGHT_MAX as i64,
        false => -i64::MAX / 10,
    };
    let deflectioncontrib = match deflection < DEFLECTION_MAX && deflection > -DEFLECTION_MAX {
        true => 0 as i64,
        false => -i64::MAX / 10,
    };
    let flight_hours_contrib =
        match flight_hours > FATIGUE_FLIGHT_HOURS_MIN && flight_hours < FATIGUE_FLIGHT_HOURS_MAX {
            true => 0 as i64,
            false => -i64::MAX / 10,
        };

    return lengthcontrib
        + heightcontrib
        + widthcontrib
        + weightcontrib
        + costcontrib
        + factorcontrib
        + deflectioncontrib
        + flight_hours_contrib;
}

fn get_rbeam_gscore(beam: beams::RectBeam) -> i64 {
    return match beam.Thickness >= 0.1 * beam.Width
        && beam.Thickness >= 0.1 * beam.Height
        && beam.Thickness <= 0.5 * beam.Height
        && beam.Thickness <= 0.5 * beam.Width
    {
        true => 0,
        false => -i64::MAX / 10,
    };
}

fn get_tbeam_gscore(beam: TBeam) -> i64 {
    return match beam.StemThickness >= 0.1 * beam.Width
        && beam.FlangeThickness >= 0.1 * beam.Height
        && beam.StemThickness <= 0.5 * beam.Height
        && beam.FlangeThickness <= 0.5 * beam.Width
    {
        true => 0,
        false => -i64::MAX / 10,
    };
}

fn get_ibeam_gscore(beam: IBeam) -> i64 {
    return match beam.CenterThickness >= 0.1 * beam.Width
        && beam.FlangeThickness >= 0.1 * beam.Height
        && beam.CenterThickness <= beam.Width
        && beam.FlangeThickness <= 0.5 * beam.Height
    {
        true => 0,
        false => -i64::MAX / 10,
    };
}

fn get_rbeam_score(beam: RectBeam) -> i64 {
    let height = beam.Height;
    let width = beam.Width;
    let length = beam.Length;

    let weight = beam.weight();
    let cost = beam.cost();

    let factor_of_safety = beam.factor_of_safety();
    let deflection = beam.vertical_deflection();
    let flight_hours = beam.flight_hours();

    return score(
        length,
        height,
        width,
        weight,
        cost,
        factor_of_safety,
        deflection,
        flight_hours,
    ) + get_rbeam_gscore(beam);
}

fn get_tbeam_score(beam: TBeam) -> i64 {
    let height = beam.Height;
    let width = beam.Width;
    let length = beam.Length;

    let weight = beam.weight();
    let cost = beam.cost();

    let factor_of_safety = beam.factor_of_safety();
    let deflection = beam.vertical_deflection();
    let flight_hours = beam.flight_hours();

    return score(
        length,
        height,
        width,
        weight,
        cost,
        factor_of_safety,
        deflection,
        flight_hours,
    ) + get_tbeam_gscore(beam);
}

fn get_ibeam_score(beam: IBeam) -> i64 {
    let height = beam.Height;
    let width = beam.Width;
    let length = beam.Length;

    let weight = beam.weight();
    let cost = beam.cost();

    let factor_of_safety = beam.factor_of_safety();
    let deflection = beam.vertical_deflection();
    let flight_hours = beam.flight_hours();

    return score(
        length,
        height,
        width,
        weight,
        cost,
        factor_of_safety,
        deflection,
        flight_hours,
    ) + get_ibeam_gscore(beam);
}

//Handling 'Fitness' for each beam.
#[derive(Eq, PartialEq, PartialOrd, Ord)]
struct BeamFitness {
    score: i64,
}

impl Fitness for BeamFitness {
    fn zero() -> BeamFitness {
        BeamFitness { score: 0 }
    }

    fn abs_diff(&self, other: &BeamFitness) -> BeamFitness {
        BeamFitness {
            score: self.score - other.score,
        }
    }
}
impl Phenotype<i64> for RectBeam {
    fn fitness(&self) -> i64 {
        return get_rbeam_score(*self);
    }

    fn crossover(&self, other: &RectBeam) -> RectBeam {
        return beams::RectBeam {
            Material: self.Material,
            Length: LENGTH,
            Height: HEIGHT,
            Width: (self.Width + other.Width) / 2.0,
            Thickness: (self.Thickness + other.Thickness) / 2.0,
        };
    }

    fn mutate(&self) -> RectBeam {
        let mut rng = thread_rng();
        let width_amp = (WIDTH_MAX - WIDTH_MIN) / 2.0;
        let thickness_amp = width_amp * 0.1;
        let new_width: f64 = self.Width + (rng.gen::<f64>() - 0.5) * width_amp / 10.0;
        let new_thickness: f64 = self.Thickness + (rng.gen::<f64>() - 0.5) * thickness_amp / 10.0;
        return beams::RectBeam {
            Material: self.Material,
            Length: LENGTH,
            Height: HEIGHT,
            Width: new_width,
            Thickness: new_thickness,
        };
    }
}

impl Phenotype<i64> for TBeam {
    fn fitness(&self) -> i64 {
        return get_tbeam_score(*self);
    }

    fn crossover(&self, other: &TBeam) -> TBeam {
        return beams::TBeam {
            Material: self.Material,
            Length: LENGTH,
            Height: HEIGHT,
            Width: (self.Width + other.Width) / 2.0,
            StemThickness: (self.StemThickness + other.StemThickness) / 2.0,
            FlangeThickness: (self.FlangeThickness + other.FlangeThickness) / 2.0,
        };
    }

    fn mutate(&self) -> TBeam {
        let mut rng = thread_rng();
        let width_amp = (WIDTH_MAX - WIDTH_MIN) / 2.0;
        let thickness_amp = width_amp * 0.1;
        let new_width: f64 = self.Width + (rng.gen::<f64>() - 0.5) * width_amp / 20.0;
        let new_stemthickness: f64 =
            self.StemThickness + (rng.gen::<f64>() - 0.5) * thickness_amp / 20.0;
        let new_flangethickness: f64 =
            self.FlangeThickness + (rng.gen::<f64>() - 0.5) * thickness_amp / 20.0;
        return beams::TBeam {
            Material: self.Material,
            Length: LENGTH,
            Height: HEIGHT,
            Width: new_width,
            StemThickness: new_stemthickness,
            FlangeThickness: new_flangethickness,
        };
    }
}

impl Phenotype<i64> for IBeam {
    fn fitness(&self) -> i64 {
        return get_ibeam_score(*self);
    }

    fn crossover(&self, other: &IBeam) -> IBeam {
        return beams::IBeam {
            Material: self.Material,
            Length: LENGTH,
            Height: HEIGHT,
            Width: (self.Width + other.Width) / 2.0,
            CenterThickness: (self.CenterThickness + other.CenterThickness) / 2.0,
            FlangeThickness: (self.FlangeThickness + other.FlangeThickness) / 2.0,
        };
    }

    fn mutate(&self) -> IBeam {
        let mut rng = thread_rng();
        let width_amp = (WIDTH_MAX - WIDTH_MIN) / 2.0;
        let thickness_amp = width_amp * 0.1;
        let new_width: f64 = self.Width + (rng.gen::<f64>() - 0.5) * width_amp / 20.0;
        let new_stemthickness: f64 =
            self.CenterThickness + (rng.gen::<f64>() - 0.5) * thickness_amp / 20.0;
        let new_flangethickness: f64 =
            self.FlangeThickness + (rng.gen::<f64>() - 0.5) * thickness_amp / 20.0;
        return beams::IBeam {
            Material: self.Material,
            Length: LENGTH,
            Height: HEIGHT,
            Width: new_width,
            CenterThickness: new_stemthickness,
            FlangeThickness: new_flangethickness,
        };
    }
}

fn get_material_string(material: Material) -> String {
    let name = match material {
        Material::Steel1018 => "Steel 1018",
        Material::StainlessSteel174PH => "Stainless Steel 17-PH",
        Material::SteelSae4340 => "Steel SAE 4340",
        Material::Aluminum7075T6 => "Aluminum 7075-T6",
        Material::Aluminum2024T4 => "Aluminum 2024-T4",
        Material::Aluminum6061T6 => "Aluminum 6061-T6",
        Material::TitaniumAlloyTi6AL4V => "Titanum Alloy Ti-6Al-4V",
    };

    return String::from(name);
}

fn get_rbeam_pop(size: u32, material: beams::Material) -> Vec<RectBeam> {
    let mut rng = thread_rng();
    return (0..size)
        .map(|i| beams::RectBeam {
            Material: material,
            Length: LENGTH,
            Height: HEIGHT,
            Width: rng.gen::<f64>() * (WIDTH_MAX - WIDTH_MIN) + WIDTH_MIN,
            Thickness: (rng.gen::<f64>() * (WIDTH_MAX - WIDTH_MIN) + WIDTH_MIN) / 10.0,
        })
        .collect();
}

fn get_tbeam_pop(size: u32, material: beams::Material) -> Vec<TBeam> {
    let mut rng = thread_rng();
    return (0..size)
        .map(|i| beams::TBeam {
            Material: material,
            Length: LENGTH,
            Height: HEIGHT,
            Width: rng.gen::<f64>() * (WIDTH_MAX - WIDTH_MIN) + WIDTH_MIN,
            StemThickness: (rng.gen::<f64>() * (WIDTH_MAX - WIDTH_MIN) + WIDTH_MIN) / 10.0,
            FlangeThickness: (rng.gen::<f64>() * (WIDTH_MAX - WIDTH_MIN) + HEIGHT) / 10.0,
        })
        .collect();
}

fn get_ibeam_pop(size: u32, material: beams::Material) -> Vec<IBeam> {
    let mut rng = thread_rng();
    return (0..size)
        .map(|i| beams::IBeam {
            Material: material,
            Length: LENGTH,
            Height: HEIGHT,
            Width: rng.gen::<f64>() * (WIDTH_MAX - WIDTH_MIN) + WIDTH_MIN,
            CenterThickness: (rng.gen::<f64>() * (WIDTH_MAX - WIDTH_MIN) + WIDTH_MIN) / 10.0,
            FlangeThickness: (rng.gen::<f64>() * (WIDTH_MAX - WIDTH_MIN) + HEIGHT) / 10.0,
        })
        .collect();
}

fn output_beam_specs(beam: impl Stress + Cost + Weight) {
    println!("Specs \n Cost: {} (<500000) \n Weight: {} (<78000) \n Flight Hours: {} (42000<->500000)\n Deflection: {} (-70<->70) \n FOS: {} (1.4<->2.2)", beam.cost(),beam.weight(),beam.flight_hours(),beam.vertical_deflection(),beam.factor_of_safety());
}

fn output_rbeam(beam: RectBeam) {
    println!(
        "Rectangular Beam\n Length: {} \n Height: {} \n Width: {} \n Thickness: {} \n Material: {}\n Score: {}", beam.Length, beam.Height, beam.Width, beam.Thickness, get_material_string(beam.Material),get_rbeam_score(beam)
    );
    output_beam_specs(beam);
}

fn output_tbeam(beam: TBeam) {
    println!(
        "T Beam\n Length: {} \n Height: {} \n Width: {} \n StemThickness: {} \n FlangeThickness: {} \n Material: {}\n Score: {}", beam.Length, beam.Height, beam.Width, beam.StemThickness,beam.FlangeThickness, get_material_string(beam.Material),get_tbeam_score(beam)
    );
    output_beam_specs(beam);
}

fn output_ibeam(beam: IBeam) {
    println!(
        "I Beam\n Length: {} \n Height: {} \n Width: {} \n CenterThickness: {} \n FlangeThickness: {} \n Material: {}\n Score: {}", beam.Length, beam.Height, beam.Width, beam.CenterThickness,beam.FlangeThickness, get_material_string(beam.Material),get_ibeam_score(beam)
    );

    output_beam_specs(beam);
}

#[derive(Copy, Clone)]
enum Beams {
    R(RectBeam),
    T(TBeam),
    I(IBeam),
}

fn main() {
    let (tx, rx): (Sender<Beams>, Receiver<Beams>) = mpsc::channel();

    for material in beams::Material::iter() {
        let transmitter = tx.clone();
        thread::spawn(move || {
            let mut rbeam_pop = get_rbeam_pop(POP_SIZE, material);
            let mut sim = Simulator::builder(&mut rbeam_pop)
                .set_selector(Box::new(UnstableMaximizeSelector::new(POP_SURVIVORS)))
                .set_max_iters(GENETIC_ITERS)
                .build();
            sim.run();
            let result = sim.get().unwrap();
            transmitter
                .send(Beams::R(*result))
                .expect("Fail to send result.");
            drop(transmitter);
        });
        let transmitter = tx.clone();
        thread::spawn(move || {
            let mut tbeam_pop = get_tbeam_pop(POP_SIZE, material);
            let mut sim = Simulator::builder(&mut tbeam_pop)
                .set_selector(Box::new(UnstableMaximizeSelector::new(POP_SURVIVORS)))
                .set_max_iters(GENETIC_ITERS)
                .build();
            sim.run();
            let result = sim.get().unwrap();
            transmitter
                .send(Beams::T(*result))
                .expect("Fail to send result.");
            drop(transmitter);
        });
        let transmitter = tx.clone();
        thread::spawn(move || {
            let mut ibeam_pop = get_ibeam_pop(POP_SIZE, material);
            let mut sim = Simulator::builder(&mut ibeam_pop)
                .set_selector(Box::new(UnstableMaximizeSelector::new(POP_SURVIVORS)))
                .set_max_iters(GENETIC_ITERS)
                .build();
            sim.run();
            let result = sim.get().unwrap();
            transmitter
                .send(Beams::I(*result))
                .expect("Fail to send result.");
            drop(transmitter);
        });
    }

    drop(tx);
    let mut best_beams: Vec<Beams> = Vec::new();

    for beam in rx {
        best_beams.push(beam);
    }

    let mut current_max_score = -i64::MAX;
    let mut best_beam = Beams::R(beams::RectBeam {
        Material: beams::Material::Steel1018,
        Height: 10.0,
        Width: 10.0,
        Thickness: 10.0,
        Length: 10.0,
    });

    for beam in best_beams {
        let score = match beam {
            Beams::R(a) => get_rbeam_score(a),
            Beams::T(a) => get_tbeam_score(a),
            Beams::I(a) => get_ibeam_score(a),
        };
        if score > current_max_score {
            current_max_score = score;
            best_beam = beam;
        }
    }

    match best_beam {
        Beams::R(a) => output_rbeam(a),
        Beams::T(a) => output_tbeam(a),
        Beams::I(a) => output_ibeam(a),
    }
}
