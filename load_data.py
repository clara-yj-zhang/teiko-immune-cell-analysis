import csv
import os
import sqlite3

DB_NAME = "immune_cells.db"
DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), DB_NAME)
CSV_NAME = "cell-count.csv"
CSV_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), CSV_NAME)

SCHEMA = """
CREATE TABLE IF NOT EXISTS subjects (
    subject     TEXT PRIMARY KEY,
    project     TEXT NOT NULL,
    condition   TEXT NOT NULL,
    age         INTEGER NOT NULL,
    sex         TEXT NOT NULL CHECK(sex IN ('M', 'F'))
);

CREATE TABLE IF NOT EXISTS treatments (
    id          INTEGER PRIMARY KEY,
    subject     TEXT NOT NULL REFERENCES subjects(subject),
    treatment   TEXT NOT NULL,
    response    TEXT CHECK(response IN ('yes', 'no', NULL))
);

CREATE TABLE IF NOT EXISTS samples (
    sample                      TEXT PRIMARY KEY,
    treatment_id                INTEGER NOT NULL REFERENCES treatments(id),
    sample_type                 TEXT NOT NULL,
    time_from_treatment_start   INTEGER NOT NULL,
    b_cell                      INTEGER NOT NULL,
    cd8_t_cell                  INTEGER NOT NULL,
    cd4_t_cell                  INTEGER NOT NULL,
    nk_cell                     INTEGER NOT NULL,
    monocyte                    INTEGER NOT NULL
);

-- Indexes for common query filters
CREATE INDEX IF NOT EXISTS idx_subjects_condition   ON subjects(condition);
CREATE INDEX IF NOT EXISTS idx_treatments_subject   ON treatments(subject);
CREATE INDEX IF NOT EXISTS idx_treatments_treatment ON treatments(treatment);
CREATE INDEX IF NOT EXISTS idx_samples_treatment_id ON samples(treatment_id);
CREATE INDEX IF NOT EXISTS idx_samples_type_time    ON samples(sample_type, time_from_treatment_start);

-- Convenience view: flattens all three tables for easy querying
CREATE VIEW IF NOT EXISTS sample_details AS
SELECT
    s.sample,
    s.sample_type,
    s.time_from_treatment_start,
    s.b_cell,
    s.cd8_t_cell,
    s.cd4_t_cell,
    s.nk_cell,
    s.monocyte,
    t.id        AS treatment_id,
    t.treatment,
    t.response,
    sub.subject,
    sub.project,
    sub.condition,
    sub.age,
    sub.sex
FROM samples s
JOIN treatments t   ON s.treatment_id = t.id
JOIN subjects sub   ON t.subject = sub.subject;
"""


def load_data():
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.executescript(SCHEMA)

    with open(CSV_PATH, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        subjects = {}
        treatments = {}   # (subject, treatment, response) -> id
        treatment_rows = []
        sample_rows = []
        treatment_id_counter = 1

        for row in reader:
            subject = row["subject"]
            if subject not in subjects:
                subjects[subject] = (
                    subject,
                    row["project"],
                    row["condition"],
                    int(row["age"]),
                    row["sex"],
                )

            treatment_key = (subject, row["treatment"], row["response"] or None)
            if treatment_key not in treatments:
                treatments[treatment_key] = treatment_id_counter
                treatment_rows.append((
                    treatment_id_counter,
                    subject,
                    row["treatment"],
                    row["response"] if row["response"] else None,
                ))
                treatment_id_counter += 1

            sample_rows.append((
                row["sample"],
                treatments[treatment_key],
                row["sample_type"],
                int(row["time_from_treatment_start"]),
                int(row["b_cell"]),
                int(row["cd8_t_cell"]),
                int(row["cd4_t_cell"]),
                int(row["nk_cell"]),
                int(row["monocyte"]),
            ))

    conn.executemany(
        "INSERT OR REPLACE INTO subjects VALUES (?,?,?,?,?)",
        subjects.values(),
    )
    conn.executemany(
        "INSERT OR REPLACE INTO treatments VALUES (?,?,?,?)",
        treatment_rows,
    )
    conn.executemany(
        "INSERT OR REPLACE INTO samples VALUES (?,?,?,?,?,?,?,?,?)",
        sample_rows,
    )
    conn.commit()
    conn.close()

    print(f"Database created at: {DB_PATH}")
    print(f"Loaded {len(subjects)} subjects, {len(treatment_rows)} treatment enrollments, and {len(sample_rows)} samples.")


if __name__ == "__main__":
    load_data()
