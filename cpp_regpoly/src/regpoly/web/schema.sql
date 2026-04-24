-- regpoly web application schema.

CREATE TABLE IF NOT EXISTS schema_version (
    version     INTEGER NOT NULL,
    applied_at  TEXT    NOT NULL DEFAULT (datetime('now'))
);

-- ================================================================
-- PRIMITIVE GENERATORS
-- ================================================================

CREATE TABLE IF NOT EXISTS primitive_search_run (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    family            TEXT    NOT NULL,
    L                 INTEGER NOT NULL,
    k                 INTEGER NOT NULL,
    structural_params TEXT    NOT NULL,
    fixed_params      TEXT    NOT NULL,
    max_tries         INTEGER,
    max_seconds       REAL,
    status            TEXT    NOT NULL DEFAULT 'pending',
    tries_done        INTEGER NOT NULL DEFAULT 0,
    found_count       INTEGER NOT NULL DEFAULT 0,
    elapsed_seconds   REAL,
    error_message     TEXT,
    created_at        TEXT    NOT NULL DEFAULT (datetime('now')),
    started_at        TEXT,
    finished_at       TEXT,
    -- Exhaustive-search mode bookkeeping.
    search_mode       TEXT    NOT NULL DEFAULT 'random',
    enum_index        INTEGER NOT NULL DEFAULT 0,
    enum_total        TEXT,
    enum_axes         TEXT
);

CREATE INDEX IF NOT EXISTS idx_psr_family ON primitive_search_run(family);
CREATE INDEX IF NOT EXISTS idx_psr_status ON primitive_search_run(status);

CREATE TABLE IF NOT EXISTS primitive_generator (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    search_run_id     INTEGER REFERENCES primitive_search_run(id) ON DELETE SET NULL,
    family            TEXT    NOT NULL,
    L                 INTEGER NOT NULL,
    k                 INTEGER NOT NULL,
    structural_params TEXT    NOT NULL,
    search_params     TEXT    NOT NULL,
    all_params        TEXT    NOT NULL,
    found_at_try      INTEGER,
    created_at        TEXT    NOT NULL DEFAULT (datetime('now')),
    -- Analysis fields (populated asynchronously once the generator
    -- has been added; NULL until the analysis worker runs).
    char_poly         TEXT,               -- hex string of the characteristic polynomial
    hamming_weight    INTEGER,            -- popcount of char_poly
    pis_gaps          TEXT,               -- JSON array of gap_v for v=1..L
    pis_se            INTEGER,            -- sum of gaps
    pis_elapsed       REAL,               -- seconds to compute
    pis_computed_at   TEXT,               -- datetime when analysis finished
    pis_error         TEXT,               -- last error message, if any
    library_id        TEXT,               -- slug from docs/library/*.yaml
                                          -- when this row is the canonical
                                          -- published generator for that id
    UNIQUE(family, structural_params, search_params)
);

-- (idx_pg_pending_analysis is created by _migrate() after the
-- analysis columns are added, to work with older DBs.)

CREATE INDEX IF NOT EXISTS idx_pg_family ON primitive_generator(family);
CREATE INDEX IF NOT EXISTS idx_pg_search_run ON primitive_generator(search_run_id);
CREATE INDEX IF NOT EXISTS idx_pg_k ON primitive_generator(k);

-- ================================================================
-- TEMPERING SEARCH
-- ================================================================

CREATE TABLE IF NOT EXISTS tempering_search_run (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    test_type         TEXT    NOT NULL,
    test_config       TEXT    NOT NULL,
    Lmax              INTEGER NOT NULL,
    nb_tries          INTEGER NOT NULL,
    optimizer_config  TEXT,
    status            TEXT    NOT NULL DEFAULT 'pending',
    combos_total      INTEGER,
    combos_done       INTEGER NOT NULL DEFAULT 0,
    best_se           INTEGER,
    elapsed_seconds   REAL,
    error_message     TEXT,
    created_at        TEXT    NOT NULL DEFAULT (datetime('now')),
    started_at        TEXT,
    finished_at       TEXT
);

CREATE INDEX IF NOT EXISTS idx_tsr_status ON tempering_search_run(status);

CREATE TABLE IF NOT EXISTS tempering_search_component (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    search_run_id     INTEGER NOT NULL REFERENCES tempering_search_run(id) ON DELETE CASCADE,
    component_index   INTEGER NOT NULL,
    shared_with_component INTEGER,
    tempering_config  TEXT    NOT NULL,
    UNIQUE(search_run_id, component_index)
);

CREATE TABLE IF NOT EXISTS tempering_search_generator (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    component_id      INTEGER NOT NULL REFERENCES tempering_search_component(id) ON DELETE CASCADE,
    generator_id      INTEGER NOT NULL REFERENCES primitive_generator(id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_tsg_component ON tempering_search_generator(component_id);

-- ================================================================
-- TESTED GENERATORS
-- ================================================================

CREATE TABLE IF NOT EXISTS tested_generator (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    search_run_id     INTEGER REFERENCES tempering_search_run(id) ON DELETE SET NULL,
    Lmax              INTEGER NOT NULL,
    k_g               INTEGER NOT NULL,
    J                 INTEGER NOT NULL,
    created_at        TEXT    NOT NULL DEFAULT (datetime('now')),
    library_id        TEXT
);

CREATE INDEX IF NOT EXISTS idx_tg_search_run ON tested_generator(search_run_id);
CREATE INDEX IF NOT EXISTS idx_tg_k_g ON tested_generator(k_g);
-- idx_pg_library_id / idx_tg_library_id are created in database.py::_migrate
-- (after ALTER TABLE adds library_id on older DBs).

CREATE TABLE IF NOT EXISTS tested_generator_component (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    tested_gen_id     INTEGER NOT NULL REFERENCES tested_generator(id) ON DELETE CASCADE,
    component_index   INTEGER NOT NULL,
    generator_id      INTEGER REFERENCES primitive_generator(id) ON DELETE SET NULL,
    family            TEXT    NOT NULL,
    L                 INTEGER NOT NULL,
    k                 INTEGER NOT NULL,
    all_params        TEXT    NOT NULL,
    tempering_params  TEXT    NOT NULL,
    UNIQUE(tested_gen_id, component_index)
);

CREATE TABLE IF NOT EXISTS test_result (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    tested_gen_id     INTEGER NOT NULL REFERENCES tested_generator(id) ON DELETE CASCADE,
    test_type         TEXT    NOT NULL,
    test_config       TEXT    NOT NULL,
    se                INTEGER,
    is_me             INTEGER,
    secf              INTEGER,
    is_cf             INTEGER,
    score             REAL,
    detail            TEXT    NOT NULL,
    elapsed_seconds   REAL,
    created_at        TEXT    NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX IF NOT EXISTS idx_tr_tested_gen ON test_result(tested_gen_id);
CREATE INDEX IF NOT EXISTS idx_tr_test_type ON test_result(test_type);
CREATE INDEX IF NOT EXISTS idx_tr_se ON test_result(se);

-- ================================================================
-- SEARCH PROGRESS (live updates)
-- ================================================================

CREATE TABLE IF NOT EXISTS search_progress (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    search_type       TEXT    NOT NULL,
    search_run_id     INTEGER NOT NULL,
    tries_done        INTEGER NOT NULL DEFAULT 0,
    found_count       INTEGER NOT NULL DEFAULT 0,
    current_info      TEXT,
    message           TEXT,
    updated_at        TEXT    NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX IF NOT EXISTS idx_sp_lookup ON search_progress(search_type, search_run_id);

-- ================================================================
-- YAML IMPORT TRACKING
-- ================================================================

CREATE TABLE IF NOT EXISTS yaml_import (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    file_path         TEXT    NOT NULL UNIQUE,
    import_type       TEXT    NOT NULL,
    row_count         INTEGER NOT NULL,
    imported_at       TEXT    NOT NULL DEFAULT (datetime('now'))
);
