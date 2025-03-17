use core::error::Error;
use core::f64;
use num::complex::Complex64;
use plotters::chart::ChartBuilder;
use plotters::prelude::{BitMapBackend, IntoDrawingArea, IntoSegmentedCoord};
use plotters::series::Histogram;
use plotters::style::{Color, RED, WHITE};
use std::{f64::consts::PI, vec};

fn dft(input: &[f64]) -> Vec<Complex64> {
    let n = input.len();
    println!("{:?}", n);

    let mut out: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); n];

    for k in 0..n {
        let mut sum = Complex64::new(0.0, 0.0);
        for l in 0..n {
            let angle: f64 = (-2.0) * PI * (k as f64) * (l as f64) / (n as f64);
            let wn = Complex64::new(
                if angle.cos().abs() < 1e-10 {
                    0.0
                } else {
                    angle.cos()
                },
                if angle.sin().abs() < 1e-10 {
                    0.0
                } else {
                    angle.sin()
                },
            );
            sum += input[l] * wn;
            sum = Complex64::new(
                if sum.re.abs() > 1e-10 { sum.re } else { 0.0 },
                if sum.im.abs() > 1e-10 { sum.im } else { 0.0 },
            );
        }
        out[k] = sum;
    }
    out
}

fn idft(spectrum: &[Complex64]) -> Vec<Complex64> {
    let n = spectrum.len();
    let mut out: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); n];

    for l in 0..n {
        let mut sum = Complex64::new(0.0, 0.0);
        for k in 0..n {
            let angle: f64 = 2.0 * PI * (k as f64) * (l as f64) / (n as f64);
            let wn = Complex64::new(angle.cos(), angle.sin());
            sum += spectrum[k] * wn;
        }
        out[l] = sum / (n as f64);
    }
    out
}

fn fft(signal: &mut [Complex64]) {
    let n = signal.len();

    if n <= 0 {
        return;
    }

    let (mut even, mut odd) = (
        vec![Complex64::new(0.0, 0.0); n / 2],
        vec![Complex64::new(0.0, 0.0); n / 2],
    );
    for i in 0..(n / 2) {
        even[i] = signal[2 * i];
        odd[i] = signal[2 * i + 1];
    }

    fft(&mut even);
    fft(&mut odd);

    for k in 0..((n as usize) / 2) {
        let mut t = Complex64::from_polar(1.0, (-2.0) * PI * (k as f64) / (n as f64)) * odd[k];
        t = Complex64::new(
            if t.re.abs() > 1e-10 { t.re } else { 0.0 },
            if t.im.abs() > 1e-10 { t.im } else { 0.0 },
        );
        signal[k] = even[k] + t;
        signal[k + n / 2] = even[k] - t;
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let num_samples = 100;

    let input = (0..num_samples)
        .map(|i| (2.0 * PI * i as f64 / num_samples as f64).sin())
        .collect::<Vec<_>>();

    println!("{:?}", dft(&input));

    let use_fft = true;

    let complex = if use_fft {
        let mut c: Vec<Complex64> = input.iter().map(|x| Complex64::new(*x, 0.0)).collect();
        fft(&mut c);
        c
    } else {
        dft(&*input)
    };

    println!("{:?}", complex);

    let root = BitMapBackend::new("./plot.png", (640 * 2, 480 * 2)).into_drawing_area();

    root.fill(&WHITE)?;

    let data = complex
        .iter()
        .map(|complex| (complex.re * complex.re + complex.im * complex.im).sqrt())
        .collect::<Vec<_>>();

    let y_max = data.iter().cloned().max_by(|a, b| a.total_cmp(b));

    let mut chart = ChartBuilder::on(&root)
        .x_label_area_size(35)
        .y_label_area_size(40)
        .margin(5)
        .caption(if use_fft { "FFT" } else { "DFT" }, ("sans-serif", 50.0))
        .build_cartesian_2d(
            (0u32..complex.len() as u32).into_segmented(),
            0f64..y_max.unwrap_or(10.0),
        )?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .bold_line_style(WHITE.mix(0.3))
        .y_desc("|Cn|")
        .x_desc("n")
        .axis_desc_style(("sans-serif", 15))
        .draw()?;

    chart.draw_series(
        Histogram::vertical(&chart)
            .style(RED.mix(0.5).filled())
            .data(data.iter().enumerate().map(|(i, x)| (i as u32, *x))),
    )?;

    Ok(())
}
