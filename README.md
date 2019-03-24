# Tatu Software

Geophysics Electromagnetic Modeling in 1D Layered Media (<a target="blank" href="https://tatusoftware.com">https://tatusoftware.com</a>)

--------------------------------------------------------------------------------

## Dependencies

- gfortran 8 or later (Fortran 2008)
- GNU Make

## How to build

Just run `make` in the Tatu root directory to create the `tatu` executable in the same directory.

The `Makefile` has two targets:

- development (default)
- production

The `development target` use gfortran flags useful to debug development process and also compile faster. To compile this target use

```shell
$ make development
```

or just

```shell
$ make
```

The `production target` use gfortran flags that will optimize the executable, making it run faster, but compile slower. To compile this target use

```shell
$ make production
```

## How to use

To simulate your models you'll need an input file in <a target="blank" href="https://json.org">JSON (JavaScript Object Notation)</a> format. The code below shows an example input file (`input.json`) with all required keys:

```json
{
  "transmitter": {
    "model": "vmd",
    "direction": "x",
    "initial": [0, 0, 0],
    "step": 0,
    "final": 0
  },
  "receiver": {
    "direction": "x",
    "initial": [100, 0, 0],
    "step": 0,
    "final": 0
  },
  "frequency": {
    "initial": 1e-1,
    "samples": 1,
    "final": 1e-1
  },
  "layers": {
    "number": 3,
    "resistivity": [100, 500, 10],
    "thickness": [1e2, 50]
  }
}
```

To simulate this model run

```shell
$ ./tatu --input-file input.json --output-file output-name
```

to get the output in JSON (JavaScript Object Notation) format in a file named `output-name.json`.

Or run

```shell
$ ./tatu --input-file input.json --output-file output-name --output-type ssv
```

to get the output in SSV (Space-Separated Values) format in a file named `output-name.ssv`.

## Input File Structure

The input file has four sections:

- transmitter
- receiver
- frequency
- layers

Theses sections will contain the information needed to simulate the model. The following shows each one's structure.

### transmitter

key       | type        | description                             | unit
--------- | ----------- | --------------------------------------- | -----
model     | string      | the source model name                   | -
direction | string      | the direction to move it at (x, y or z) | -
initial   | real array  | initial position 3D coordinates         | meter
step      | real number | position increment step at direction    | meter
final     | real number | final position at direction             | meter

At the moment, `tatu` supports the following transmitter models:

- hedx (horizontal electric dipole in x direction)
- hedy (horizontal electric dipole in y direction)
- ved (vertical electric dipole)
- hmdx (horizontal magnetic dipole in x direction)
- hmdy (horizontal magnetic dipole in y direction)
- vmd (vertical magnetic dipole)

### receiver

key       | type        | description                             | unit
--------- | ----------- | --------------------------------------- | -----
direction | string      | the direction to move it at (x, y or z) | -
initial   | real array  | initial position 3D coordinates         | meter
step      | real number | position increment step at direction    | meter
final     | real number | final position at direction             | meter

### frequency

key     | type           | description    | unit
------- | -------------- | -------------- | -----
initial | real number    | initial value  | Hertz
samples | integer number | samples number | -
final   | real number    | final value    | Hertz

The frequency samples are logarithmically interpolated.

### layers

key         | type           | description        | unit
----------- | -------------- | ------------------ | -----------
number      | integer number | layers number      | -
resistivity | real array     | layers resistivity | Ohm x meter
thickness   | real array     | layers thickness   | meter

The resistivity array length has to be equal the layers number.

The thickness array length has to be equal the layers number minus one.

## Output Files Structure

At the moment, the `tatu` software can save output data in two different formats:

- JSON (JavaScript Object Notation)
- SSV (Space-Separated Values)

### JSON format

Output in JSON format contains three sections:

- output (see below)
- input (the input model)
- unique (transmitter, receiver and frequency unique values)

The output section has two subsections:

- labels (array with ordered values labels)
- values (array of values arrays)

The labels array informs the order of the values in each item in values array.

```json
{
  "output": {
    "labels": [
      "transmitter",
      "frequency",
      "receiver",
      "ExReal",
      "ExImag",
      "EyReal",
      "EyImag",
      "EzReal",
      "EzImag",
      "HxReal",
      "HxImag",
      "HyReal",
      "HyImag",
      "HzReal",
      "HzImag"
    ],
    "values": [
      [
        0.0E+0,
        0.10000000000000001E+0,
        0.1E+3,
        0.0E+0,
        0.0E+0,
        -0.28521355723864792E-15,
        -0.6283174409521896E-11,
        0.69358017024460419E+0,
        0.6952704522069525E+0,
        -0.16264080201972775E-12,
        0.22258411762892078E-11,
        0.0E+0,
        0.0E+0,
        -0.79577978723228242E-7,
        -0.5556887181527204E-11
      ]
    ]
  },
  "input": {
    "transmitter": {
      "model": "vmd  ",
      "direction": "x",
      "initial": [ 0.0E+0, 0.0E+0, 0.0E+0 ],
      "step": 0.0E+0,
      "final": 0.0E+0
    },
    "receiver": {
      "direction": "x",
      "initial": [ 0.1E+3, 0.0E+0, 0.0E+0 ],
      "step": 0.0E+0,
      "final": 0.0E+0
    },
    "frequency": {
      "initial": 0.10000000000000001E+0,
      "samples": 1,
      "final": 0.10000000000000001E+0
    },
    "layers": {
      "number": 3,
      "resistivity": [ 0.1E+3, 0.5E+3, 0.1E+2 ],
      "thickness": [ 0.1E+3, 0.5E+2 ]
    }
  },
  "unique": {
    "transmitter": [
      [ 0.0E+0, 0.0E+0, 0.0E+0 ]
    ],
    "frequency": [
      0.10000000000000001E+0
    ],
    "receiver": [
      [ 0.1E+3, 0.0E+0, 0.0E+0 ]
    ]
  }
}
```

## SSV format

Output in SSV format is table-like, in other words presents data in rows and columns. The data in each row follows the same order as shown in JSON format's labels array.

```
0.00000000000000 0.100000000000000 100.000000000000 0.00000000000000 0.00000000000000 -0.285213557238648E-015 -0.628317440952190E-011   0.00000000000000 0.00000000000000 -0.162640802019728E-012 0.222584117628921E-011 0.00000000000000 0.00000000000000 -0.795779787232282E-007 -0.555688718152720E-011
```

## Options available

To see all options available use

```shell
$ ./tatu --help
```

The output will be something like:

```shell
Tatu - Geophysics Electromagnetic Modeling in 1D Layered Media

Usage: tatu [options]

Options (*required):

    -v, --version                   Show version and exit
    -h, --help                      Show this help message
  * -i, --input-file <FILEPATH>     File to read the input data
  * -o, --output-file <FILEPATH>    File to write the output data
    -t, --output-type <FILETYPE>    Output file type: json (default), ssv or all
```
