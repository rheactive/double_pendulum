use macroquad::prelude::*;

use std::fs;
use std::io::Write;
use std::path::PathBuf;

use std::collections::LinkedList;
use std::thread::sleep;
use std::time::Duration;

use std::vec::Vec;

const FT_DESIRED: f64 = 0.01666666666667;
const PI: f64 = std::f64::consts::PI;
const SPEED: f64 = 1.0;
// trail parameter
const MEM: usize = 50;
// g
const G: f64 = 9.81;

struct MyVec {
    vec: Vec<f64>,
    n: usize,
}

impl MyVec {
    fn _new_emp(n: usize) -> MyVec {
        MyVec {
            vec: vec![0.0; n],
            n,
        }
    }

    fn new_from(a: &Vec<f64>) -> MyVec {
        MyVec {
            vec: a.to_vec(),
            n: a.len(),
        }
    }

    fn add(&self, a: &MyVec) -> MyVec {
        if self.n != a.n {
            panic!("Vector addition panic: lengths not equal")
        }
        let mut c = vec![0.0; self.n];
        for j in 0..self.n {
            c[j] = self.vec[j] + a.vec[j]
        }
        MyVec::new_from(&c)
    }

    fn scale(&self, s: f64) -> MyVec {
        let mut c = vec![0.0; self.n];
        for j in 0..self.n {
            c[j] = self.vec[j] * s
        }
        MyVec::new_from(&c)
    }
}

fn runge_kutta(
    vars: &MyVec,
    pars: &Vec<f64>,
    rhs: &dyn Fn(&MyVec, &Vec<f64>) -> MyVec,
    dt: f64,
) -> MyVec {
    let rk_1 = rhs(vars, pars);
    let rk_2 = rhs(&vars.add(&rk_1.scale(dt / 2.0)), pars);
    let rk_3 = rhs(&vars.add(&rk_2.scale(dt / 2.0)), pars);
    let rk_4 = rhs(&vars.add(&rk_3.scale(dt)), pars);

    let vars_new = vars
        .add(&rk_1.scale(dt / 6.0))
        .add(&rk_2.scale(dt / 3.0))
        .add(&rk_3.scale(dt / 3.0))
        .add(&rk_4.scale(dt / 6.0));
    vars_new
}

// pendulum struct

#[derive(Clone)]
struct Pendulum {
    period: f64,
    omega2: f64,
    angle_deg: f64,
    angle: f64,
    momentum: f64,
    angle_dot: f64,
    momentum_dot: f64,
    length: f64,
    base_vec: [f64; 2],
    ball_vec: [f64; 2],
    friction: f64,
}

impl Pendulum {
    fn new(base_vec: [f64; 2], angle_deg: f64, period: Option<f64>, friction: Option<f64>) -> Pendulum {
        let period = period.unwrap_or(2.0);
        let omega = 2.0 * PI / period;
        let omega2 = omega.powi(2);
        let angle = angle_deg / 180.0 * PI;
        let length = G / omega2;
        Pendulum {
            period,
            omega2,
            angle_deg,
            angle,
            momentum: 0.0,
            angle_dot: 0.0,
            momentum_dot: 0.0,
            length,
            base_vec,
            ball_vec: [
                base_vec[0] + length * angle.sin(),
                base_vec[1] - length * angle.cos(),
            ],
            friction: friction.unwrap_or(0.0),
        }
    }
}

// get a vector of variables from the two pendula

struct DynModel {
    pendula: Vec<Pendulum>,
    vars: MyVec,
    pars: Vec<f64>,
}

impl DynModel {
    fn new_dp (top_angle_deg: f64, bottom_angle_deg: f64) -> DynModel {
        let mut pendula = Vec::new();
        pendula.push(Pendulum::new([0.0; 2], top_angle_deg, None, None));
        pendula.push(Pendulum::new(pendula[0].ball_vec, bottom_angle_deg, None, None));
        let vars = DynModel::get_pen_vars(&pendula);
        let pars = DynModel::get_pen_pars(&pendula);

        DynModel {
            pendula,
            vars,
            pars,
        }
    }

    fn set_momenta_dp(&mut self) {
        let cs = (self.pendula[0].angle - self.pendula[1].angle).cos();
        self.pendula[0].momentum = 2.0 * self.pendula[0].angle_dot + self.pendula[1].angle_dot * cs;
        self.pendula[1].momentum = self.pendula[1].angle_dot + self.pendula[0].angle_dot * cs;
        self.vars = DynModel::get_pen_vars(&self.pendula);
        self.pars = DynModel::get_pen_pars(&self.pendula);
    }

    fn update_dp(&mut self, dt: f64) {
        self.vars = runge_kutta(&self.vars, &self.pars, &rhs_dp, dt);
    }

    fn get_pen_vars(pendula: &Vec<Pendulum>) -> MyVec {
        let mut vars = Vec::new();
        for pend in pendula {
            vars.push(pend.angle)
        }
        for pend in pendula {
            vars.push(pend.momentum)
        }
        MyVec::new_from(&vars)
    }
    
    fn get_pen_pars(pendula: &Vec<Pendulum>) -> Vec<f64> {
        let mut pars = Vec::new();
        for pend in pendula {
            pars.push(pend.omega2)
        }
        for pend in pendula {
            pars.push(pend.friction)
        }
        pars
    }
}



fn rhs_dp(vars: &MyVec, pars: &Vec<f64>) -> MyVec {
    let theta = vars.vec[0];
    let phi = vars.vec[1];
    let p = vars.vec[2];
    let q = vars.vec[3];
    let omega21 = pars[0];
    let omega22 = pars[1];
    let c_frict1 = pars[2];
    let c_frict2 = pars[3];

    let cs = (theta - phi).cos();
    let sn = (theta - phi).sin();
    let a = 1.0 / (1.0 + sn * sn);
    let theta_dot = a * (p - q * cs);
    let phi_dot = a * (2.0 * q - p * cs);
    let b = theta_dot * phi_dot * sn;
    let p_dot = -b - 2.0 * omega21 * theta.sin() - c_frict1 * p;
    let q_dot = b - omega22 * phi.sin() - c_frict2 * q;

    let v = vec![theta_dot, phi_dot, p_dot, q_dot];
    MyVec::new_from(&v)

}

fn total_energy_dp(vars: &MyVec, pars: &Vec<f64>) -> f64 {
    let theta = vars.vec[0];
    let phi = vars.vec[1];
    let p = vars.vec[2];
    let q = vars.vec[3];
    let omega21 = pars[0];
    let omega22 = pars[1];

    let cs = (theta - phi).cos();
    let sn = (theta - phi).sin();
    let a = 1.0 / (1.0 + sn * sn);
    let theta_dot = a * (p - q * cs);
    let phi_dot = a * (2.0 * q - p * cs);

    let h = 0.5 * (p * theta_dot + q * phi_dot) - (2.0 * omega21 * theta.cos() + omega22 * phi.cos());
    h
}

fn get_window_par(height: f32, width: f32) -> [f32; 4] {
    let l = height / 5.0;
    let w = width / 100.0;
    let x0 = width / 2.0;
    let y0 = height / 2.0;
    [l, w, x0, y0]
}

fn find_pendulum(theta: f64, phi: f64, par: [f32; 4]) -> [f32; 4] {
    // calculate pendulum positions
    let l = par[0];
    let x0 = par[2];
    let y0 = par[3];

    let x1 = x0 + l * (theta.sin() as f32);
    let y1 = y0 + l * (theta.cos() as f32);
    let x2 = x1 + l * (phi.sin() as f32);
    let y2 = y1 + l * (phi.cos() as f32);

    [x1, y1, x2, y2]
}

fn draw_pendulum(coords: [f32; 4], par: [f32; 4]) {
    // plot pendulum positions
    let w = par[1];
    let x0 = par[2];
    let y0 = par[3];

    let x1 = coords[0];
    let y1 = coords[1];
    let x2 = coords[2];
    let y2 = coords[3];

    draw_line(x0, y0, x1, y1, w / 2.0, BLACK);
    draw_line(x1, y1, x2, y2, w / 2.0, BLACK);
    draw_circle(x0, y0, 0.5 * w, BLACK);
    draw_circle(x1, y1, w, BLACK);
    draw_circle(x1, y1, 0.8 * w, DARKPURPLE);
    draw_circle(x2, y2, w, BLACK);
    draw_circle(x2, y2, 0.8 * w, LIGHTGRAY);
}

fn draw_trail(xy_mem: &LinkedList<(f32, f32)>, xy_last: (f32, f32), w: f32) {
    let mut iter1 = xy_mem.iter();
    let mut iter2 = xy_mem.iter();
    iter2.next();
    // calculate and plot pendulum positions
    for _jm in 0..(MEM - 1) {
        let xy1 = iter1.next().unwrap();
        let xy2 = iter2.next().unwrap();
        //draw_circle(xy1.0, xy1.1, w / 5.0, MAROON);
        draw_line(xy1.0, xy1.1, xy2.0, xy2.1, w / 5.0, DARKPURPLE);
    }
    draw_circle(xy_last.0, xy_last.1, w / 3.0, DARKPURPLE);
}

#[macroquad::main("Double Pendulum")]
async fn main() {
    //make directory
    //let dir_name = "dp_results";
    //fs::create_dir_all(dir_name).expect("Error creating directory");

    // file for data
    let fl_name = "double_pendulum.dat";
    let file_path: PathBuf = [fl_name].iter().collect();
    let mut my_file = fs::File::create(file_path).expect("Error creating file");

    // column names
    writeln!(my_file, "t theta phi p q H").expect("Error writing to file");

    loop {
        let mut top_angle_deg = 90.0;
        let mut bottom_angle_deg = 180.0;
        let mut model = DynModel::new_dp(top_angle_deg, bottom_angle_deg);

        // time parameters
        let tf = 1000.0;
        let dt_millis: u64 = 3;
        let dt = 0.001 * (dt_millis as f64);
        // frame time
        let mut ft: f64;
        let mut time_start: f64 = 0.0;
        let mut time_end: f64;

        // initial conditions
        let mut t_0 = 0.0;
        let mut t_draw = 0.0;

        //total energy
        let mut h = total_energy_dp(&model.vars, &model.pars);

        // draw the pendulum initial state

        while !is_key_down(KeyCode::Space) {
            clear_background(SKYBLUE);

            let par = get_window_par(screen_height(), screen_width());
            let coords = find_pendulum(model.vars.vec[0], model.vars.vec[1], par);

            draw_pendulum(coords, par);

            let w = par[1];

            draw_text(
                format!("Use left/right or up/down arrows to move the pendulum.").as_str(),
                5.0 * w,
                3.0 * w,
                2.5 * w,
                BLACK,
            );

            draw_text(
                format!("Angle 1: {:.1} degrees.", angle_in_degrees(top_angle_deg)).as_str(),
                5.0 * w,
                6.0 * w,
                2.5 * w,
                BLACK,
            );

            draw_text(
                format!("Angle 2: {:.1} degrees.", angle_in_degrees(bottom_angle_deg)).as_str(),
                5.0 * w,
                9.0 * w,
                2.5 * w,
                BLACK,
            );

            draw_text(
                format!("Press [F]/[A] to increase/decrease friction.").as_str(),
                5.0 * w,
                2.0 * par[3] - 9.0 * w,
                2.5 * w,
                BLACK,
            );

            draw_text(
                format!("Press [SPACE] to start.").as_str(),
                5.0 * w,
                2.0 * par[3] - 3.0 * w,
                2.5 * w,
                BLACK,
            );

            ft = get_frame_time() as f64;

            if is_key_down(KeyCode::Left) {
                top_angle_deg = top_angle_deg - SPEED * ft;
            }

            if is_key_down(KeyCode::Right) {
                top_angle_deg = top_angle_deg + SPEED * ft;
            }

            if is_key_down(KeyCode::Down) {
                bottom_angle_deg = bottom_angle_deg - SPEED * ft;
            }

            if is_key_down(KeyCode::Up) {
                bottom_angle_deg = bottom_angle_deg + SPEED * ft;
            }

            next_frame().await;
        }

        // window parameters
        let mut par: [f32; 4];
        let mut w: f32;
        let mut coords: [f32; 4];

        par = get_window_par(screen_height(), screen_width());
        coords = find_pendulum(model.vars.vec[0], model.vars.vec[1], par);

        // trail data
        let mut xy_mem = LinkedList::from([(coords[2], coords[3]); MEM]);

        //start the calculation

        while t_0 < tf && !is_key_down(KeyCode::R) {
            t_0 = t_0 + dt;

            ft = get_frame_time() as f64;

            if t_0 - t_draw > ft {
                time_end = get_time();

                if time_end - time_start < FT_DESIRED {
                    sleep(Duration::from_secs_f64(FT_DESIRED - time_end + time_start));
                }

                time_start = get_time();

                t_draw = t_0;

                clear_background(SKYBLUE);

                par = get_window_par(screen_height(), screen_width());
                w = par[1];

                coords = find_pendulum(model.vars.vec[0], model.vars.vec[1], par);

                let xy_last = xy_mem.pop_back().unwrap();

                xy_mem.push_front((coords[2], coords[3]));

                draw_trail(&xy_mem, xy_last, w);
                draw_pendulum(coords, par);

                draw_text(
                    format!("Time elapsed, seconds: {:.3}", t_0).as_str(),
                    5.0 * w,
                    6.0 * w,
                    2.0 * w,
                    BLACK,
                );
                draw_text(
                    format!("Total energy, arb. units: {:.5}", h).as_str(),
                    5.0 * w,
                    9.0 * w,
                    2.0 * w,
                    BLACK,
                );
                draw_text(
                    format!("Press [R] to reset the simulation!").as_str(),
                    par[2] + 10.0 * w,
                    6.0 * w,
                    2.0 * w,
                    BLACK,
                );

                draw_text(
                    format!("FPS is {}", get_fps()).as_str(),
                    5.0 * w,
                    2.0 * par[3] - 3.0 * w,
                    2.0 * w,
                    BLACK,
                );

                writeln!(
                    my_file,
                    "{} {} {} {} {} {}",
                    t_0,
                    angle_in_degrees(model.vars.vec[0]),
                    angle_in_degrees(model.vars.vec[1]),
                    model.vars.vec[2],
                    model.vars.vec[3],
                    h
                )
                .expect("Error writing to file");

                next_frame().await;
            }

            model.update_dp(dt);

            h = total_energy_dp(&model.vars, &model.pars)
        }
    }
}

// express angle in degrees
fn angle_in_degrees(angle: f64) -> f64 {
    let in_degrees = angle * 180.0 / PI;
    if in_degrees < 0.0 {
        360.0 + in_degrees
    } else if in_degrees > 360.0 {
        -360.0 + in_degrees
    } else {
        in_degrees
    }
}
