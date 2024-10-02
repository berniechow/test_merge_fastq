version 1.0

struct InputOutput {
  Array[String]+ reads
  String filename
}

workflow merge_fastqs {
  meta {
    version: "0.9.0"
    author: "Yunhai Luo"
    description: "Concatenate gzipped input fastqs into one set of gzipped fastq."
  }

  input {
    Array[File]+ input_files
    String sample_name
  }

  parameter_meta {
    input_files: {
      help: "[Required] An array of paths for gzipped fastqs (.fastq.gz). These fastqs will merged as a single set of R1/R2/I1/I2 when applicable. This workflow will try to sort and match R1/R2/I1/I2 in these fastqs based on their file names. Keys in file names used for sorting and matching are literal strings `_R1`, `_R2`, `_I1` or `_I2`. It will be helpful if file names follow the Illumina naming convention (https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm)."
    }
    sample_name: {
      help: "[Required] A sample name string which will be used to name the merged gzipped fastqs as `<sample_name>_(R1|R2|I1|I2)_001.fastq.gz`."
    }
  }

  call check_fastqs {
    input:
      input_file_list = write_json(input_files),
      sample_name = sample_name
  }

  scatter(io in check_fastqs.input_outputs) {
    call run_merge_fastqs {
      input:
        fastqs = io.reads,  # !FileCoercion
        merged_filename = io.filename
    }
  }

  output {
    Array[File]+ merged_fastqs = run_merge_fastqs.merged_fastq
  }
}

task check_fastqs {
  input {
    # Putting file paths in a file will help
    # 1) Avoid unnecessary localization but be careful when processing them
    #    below.
    # 2) Bypass the following croo bug where a task having one same file as
    #    both input and output will trigger a "Detected a cyclic link in DAG."
    #    error: https://github.com/ENCODE-DCC/croo/issues/40
    File input_file_list
    String sample_name 
  }

  command <<<
    python <<CODE
    import json
    import re
    from itertools import zip_longest
    from pathlib import Path

    LOWER_REGEX_READ_TYPE = {
        "(r1|_1|_r1_0\d\d).(fastq|fq).gz$": "R1",
        "(r2|_2|_r2_0\d\d).(fastq|fq).gz$": "R2",
        "(i1|_i1_0\d\d).(fastq|fq).gz$": "I1",
        "(i2|_i2_0\d\d).(fastq|fq).gz$": "I2",
    }

    # Organize input fastqs
    with open("~{input_file_list}") as fp:
        fname_sorted_input = sorted(json.load(fp), key=lambda p: Path(p).name)
    ios_dict = {
        read_type: {
            "reads": [],
            "fname_no_suffix": [],  # collect for sanity check below
            "filename": f"~{sample_name}_{read_type}.fastq.gz"
        }
        for read_type in LOWER_REGEX_READ_TYPE.values()
    }
    for f in fname_sorted_input:
        fname = Path(f).name
        for lower_regex, read_type in LOWER_REGEX_READ_TYPE.items():
            suffix_match = re.search(lower_regex, fname.lower())
            if suffix_match:
                ios_dict[read_type]["reads"].append(f)
                ios_dict[read_type]["fname_no_suffix"].append(
                    fname[:suffix_match.start()]
                )
                break
        else:
            raise ValueError(
                f"Cannot tell sequence type of {f} as "
                f"{list(LOWER_REGEX_READ_TYPE.values())}"
            )

    # Sanity checks
    filtered_ios_dict = {
        r: io
        for r, io in ios_dict.items()
        if io["reads"]
    }
    if not filtered_ios_dict:
        raise ValueError(
            f"No fastqs found, with pattern {list(LOWER_REGEX_READ_TYPE)},"
            " to be merged."
        )
    if len(filtered_ios_dict) > 1:
        # Strict check on read match: file names are expected to be otherwise
        # identical except for read type label.
        for normalized_fname_set in zip_longest(
            *[
                list(zip_longest(io["reads"], io["fname_no_suffix"]))
                for r, io in filtered_ios_dict.items()
            ]
        ):
            if None in normalized_fname_set:
                raise ValueError(
                    f"Inconsistent file count for read sets {filtered_ios_dict}"
                )
            if len({ftuple[1] for ftuple in normalized_fname_set}) != 1:
                raise ValueError(
                    "Found mismatch "
                    f"{[ftuple[0] for ftuple in normalized_fname_set]} "
                    f"in {filtered_ios_dict}"
                )

    # Serialize organized input/output
    with open("__serialized_input_output.json", "w") as fout:
        json.dump(
            [
                {
                    "reads": io["reads"],
                    "filename": io["filename"],
                }
                for io in filtered_ios_dict.values()
            ],
            fout
        )
    CODE
  >>>

  runtime {
    cpu: "1"
    memory: "2 GB"
    docker: "python:3.8.2"
  }

  output {
    Array[InputOutput]+ input_outputs = read_json("__serialized_input_output.json")
  }
}

task run_merge_fastqs {
  # Simple case: concat all input fastqs together as one in the order of input
  input {
    Array[File] fastqs
    String merged_filename
  }

  command <<<
    set -e
    set -u
    set -o pipefail
    set -x
    zcat ~{sep=" " fastqs} | gzip -n > ~{merged_filename}
  >>>

  runtime {
    cpu: "1"
    memory: "2 GB"
    docker: "python:3.8.2"
  }

  output {
    File merged_fastq = merged_filename
  }
}
