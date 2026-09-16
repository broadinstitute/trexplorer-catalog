"""Checks that the variation_clusters_and_all_repeats catalog and the legacy
variation_clusters_and_isolated_TRs catalog agree with each other.

The two files are built by separate scripts from the same variation clusters TSV and repeat
catalog, so nothing keeps them in sync except this check. The all_repeats catalog gives every
repeat its own row, plus one row per variation cluster; the legacy catalog gives a row only to
repeats that are not inside a cluster, plus the same per-cluster rows. So a repeat should have a
standalone row in the legacy catalog if and only if it has an unflagged row in the all_repeats
catalog and is not listed as a member of any variation cluster there. A repeat flagged VC_FILTER
in the all_repeats catalog is expected to be absent from the legacy catalog entirely, since the
legacy generator drops DEPTH/EXTENSION loci rather than keeping a flagged row for them
(https://github.com/PacificBiosciences/trgt-lps/issues/5).
"""

import argparse
import collections
import gzip
import re


# Loci that have a standalone row in the legacy catalog but no row at all in the all_repeats
# catalog. The two generators read different universes: the legacy generator emits an isolated
# repeat straight from the variation clusters TSV, while the all_repeats generator emits repeat
# rows only for loci that are in the repeat catalog. These loci are in the TSV but not in the
# repeat catalog, so only the legacy generator can emit them. All of them are homopolymers.
# Every id here must still be unmatched at each run: an id that stops being unmatched is a stale
# entry that would otherwise silently suppress a real inconsistency, so it fails this check.
LOCI_ABSENT_FROM_REPEAT_CATALOG = frozenset([
    "1-16213728-16213772-T",
    "1-51173963-51174007-A",
    "1-73718808-73718849-A",
    "1-73763779-73763834-A",
    "1-82429692-82429751-T",
    "1-195618934-195619004-A",
    "1-243902601-243902668-T",
    "2-11466838-11466890-C",
    "2-103340538-103340578-T",
    "2-155087927-155087993-A",
    "2-175308746-175308786-T",
    "2-237169892-237169936-A",
    "2-241409019-241409078-A",
    "3-2658950-2659010-A",
    "3-85803498-85803534-A",
    "3-107271166-107271224-A",
    "3-111355337-111355379-T",
    "3-121173323-121173350-A",
    "3-139469663-139469735-T",
    "4-3074876-3074968-C",
    "4-41745975-41746022-C",
    "4-142579857-142579937-A",
    "4-154587735-154587823-A",
    "4-157415688-157415736-T",
    "5-3424567-3424635-A",
    "5-3712600-3712668-T",
    "5-82800267-82800335-T",
    "5-146878727-146878759-G",
    "6-7692786-7692828-A",
    "6-12859987-12860026-T",
    "6-16327633-16327724-T",
    "6-24306494-24306574-A",
    "6-104226081-104226125-T",
    "6-112557951-112558018-A",
    "7-43327961-43327991-T",
    "7-144026620-144026691-T",
    "7-154776199-154776260-A",
    "7-156680376-156680424-G",
    "8-32209688-32209777-T",
    "8-143143770-143143846-A",
    "9-1651060-1651114-A",
    "9-30216224-30216263-T",
    "9-76649591-76649646-A",
    "9-129188807-129188870-T",
    "10-6679975-6680010-G",
    "10-78952323-78952359-G",
    "10-124470648-124470693-A",
    "10-127955286-127955322-A",
    "11-17064691-17064763-T",
    "11-43588594-43588637-T",
    "11-133365024-133365082-C",
    "12-5983968-5984044-A",
    "12-47390622-47390670-A",
    "12-52803344-52803372-G",
    "12-82591610-82591666-A",
    "13-54532003-54532047-T",
    "13-89042714-89042762-T",
    "13-100851724-100851772-A",
    "13-105480503-105480571-A",
    "14-23321471-23321492-G",
    "14-47612191-47612257-A",
    "14-91736479-91736529-A",
    "14-93913487-93913562-A",
    "15-36701080-36701176-A",
    "15-89042571-89042631-T",
    "16-14241963-14242011-A",
    "16-23738190-23738246-T",
    "16-82877924-82877997-A",
    "16-87604282-87604329-G",
    "17-6429800-6429868-A",
    "17-13213540-13213570-A",
    "17-48595772-48595831-T",
    "17-57222578-57222630-A",
    "18-7286408-7286460-A",
    "18-59889955-59889999-A",
    "18-71246315-71246363-T",
    "19-44827033-44827092-A",
    "20-5839486-5839526-A",
    "20-7345632-7345695-T",
    "20-9558007-9558071-T",
    "20-15317482-15317526-T",
    "20-58089575-58089659-A",
    "21-19181972-19182099-T",
    "21-20274569-20274617-G",
    "21-46617771-46617834-T",
    "X-147912049-147912111-G",
    "Y-2954699-2954743-A",
    "Y-3263110-3263158-A",
    "Y-3263416-3263486-T",
    "Y-3411660-3411696-A",
    "Y-3575723-3575759-T",
    "Y-3772409-3772421-C",
    "Y-3772789-3772837-T",
    "Y-3811618-3811658-A",
    "Y-4402918-4402978-A",
    "Y-5446867-5446906-A",
    "Y-5547245-5547289-T",
    "Y-6004215-6004241-A",
    "Y-6542908-6542965-A",
    "Y-6724574-6724613-A",
    "Y-6993189-6993257-T",
    "Y-7043527-7043571-A",
    "Y-7185317-7185385-A",
    "Y-7547583-7547623-A",
    "Y-7568215-7568255-A",
    "Y-7577817-7577853-G",
    "Y-7577969-7577997-A",
    "Y-7640314-7640375-T",
    "Y-7775467-7775577-G",
    "Y-7846943-7846983-T",
    "Y-7862389-7862473-A",
    "Y-7999838-7999902-A",
    "Y-8258258-8258303-T",
    "Y-8258316-8258341-T",
    "Y-8350055-8350139-T",
    "Y-8356114-8356158-A",
    "Y-8501335-8501359-T",
    "Y-8511087-8511163-T",
    "Y-8519498-8519594-G",
    "Y-8525492-8525528-A",
    "Y-8598153-8598197-A",
    "Y-8633403-8633423-T",
    "Y-8687938-8687978-T",
    "Y-8781943-8782055-A",
    "Y-8874697-8874748-A",
    "Y-8954510-8954553-A",
    "Y-9582529-9582583-T",
    "Y-9816671-9816778-A",
    "Y-9936779-9936873-T",
    "Y-10018088-10018118-T",
    "Y-11819430-11819490-A",
    "Y-11839880-11839968-G",
    "Y-11862987-11863031-A",
    "Y-11942439-11942531-A",
    "Y-11980424-11980763-T",
    "Y-11981553-11981613-T",
    "Y-11981653-11981689-T",
    "Y-11982088-11982132-T",
    "Y-12345803-12345839-T",
    "Y-12366980-12367012-G",
    "Y-12384502-12384538-G",
    "Y-12649169-12649229-A",
    "Y-13318158-13318224-A",
    "Y-13541541-13541565-T",
    "Y-13585211-13585243-A",
    "Y-13640727-13640835-T",
    "Y-13797170-13797214-T",
    "Y-13940665-13940733-A",
    "Y-13987207-13987253-T",
    "Y-14022415-14022455-A",
    "Y-14055475-14055521-A",
    "Y-14115982-14116050-A",
    "Y-14124353-14124385-A",
    "Y-14186619-14186663-A",
    "Y-14305460-14305500-T",
    "Y-14396603-14396627-A",
    "Y-14452470-14452528-A",
    "Y-14472030-14472070-T",
    "Y-14578222-14578261-A",
    "Y-14611543-14611591-T",
    "Y-14635531-14635551-T",
    "Y-14635562-14635602-T",
    "Y-14945046-14945082-A",
    "Y-15035513-15035568-T",
    "Y-15035625-15035675-T",
    "Y-15035685-15035697-T",
    "Y-15163066-15163162-A",
    "Y-15188020-15188088-A",
    "Y-15188101-15188119-G",
    "Y-15188126-15188138-A",
    "Y-15195628-15195684-A",
    "Y-15314131-15314186-T",
    "Y-15329869-15329953-T",
    "Y-15541824-15541848-T",
    "Y-15669603-15669694-T",
    "Y-15800434-15800470-A",
    "Y-15859770-15859802-A",
    "Y-15936843-15936876-T",
    "Y-16013850-16013931-A",
    "Y-16050586-16050626-G",
    "Y-16050661-16050695-G",
    "Y-16053587-16053629-T",
    "Y-16107056-16107080-A",
    "Y-16281075-16281155-T",
    "Y-16302501-16302577-A",
    "Y-16453354-16453416-T",
    "Y-16470035-16470118-A",
    "Y-16473881-16473956-A",
    "Y-16579095-16579135-A",
    "Y-16631634-16631756-T",
    "Y-16760350-16760423-T",
    "Y-16812722-16812758-T",
    "Y-16860294-16860326-A",
    "Y-16937590-16937656-A",
    "Y-16969637-16969673-T",
    "Y-17234447-17234459-T",
    "Y-17234482-17234494-T",
    "Y-17246969-17247009-T",
    "Y-17260392-17260448-A",
    "Y-17260490-17260514-T",
    "Y-17260549-17260561-A",
    "Y-18680512-18680687-A",
    "Y-19206306-19206371-T",
    "Y-19224281-19224311-T",
    "Y-19244116-19244209-A",
    "Y-19302872-19302944-T",
    "Y-19315131-19315143-C",
    "Y-19315187-19315196-C",
    "Y-19315280-19315387-T",
    "Y-19315879-19315981-A",
    "Y-19358337-19358389-A",
    "Y-19458591-19458746-A",
    "Y-19494950-19495000-A",
    "Y-19647923-19647959-T",
    "Y-19930715-19930763-T",
    "Y-19937747-19937795-T",
    "Y-20400677-20400713-A",
    "Y-20472970-20473014-T",
    "Y-20607708-20607741-T",
    "Y-20607772-20607784-T",
    "Y-20815639-20815699-A",
    "Y-20901812-20901852-T",
    "Y-21072785-21072889-T",
    "Y-21599708-21599794-A",
    "Y-21681708-21681748-T",
    "Y-21806014-21806055-A",
    "Y-22218922-22218996-A",
    "Y-22219030-22219078-A",
    "Y-22267099-22267123-T",
    "Y-22270954-22271019-A",
    "Y-22339545-22339610-T",
])


def open_maybe_gzipped(path):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def locus_id_of(info_field):
    return re.search(r"ID=([^;]+)", info_field).group(1)


def motifs_of(info_field):
    return re.search(r"MOTIFS=([^;]*)", info_field).group(1)


def vc_filter_of(info_field):
    match = re.search(r";VC_FILTER=([^;]*)", info_field)
    return match.group(1) if match else None


def is_variation_cluster(info_field):
    return locus_id_of(info_field).startswith("VC:")


def members_of_cluster(info_field):
    """Returns the list of member locus ids named in a variation cluster row's STRUC field."""
    return re.search(r"STRUC=<VC:([^>]*)>", info_field).group(1).split(",")


def read_catalog(path):
    """Reads a TRGT catalog BED file.

    Returns:
        tuple: (dict of cluster id -> (chrom, start, end, motifs, struc),
            dict of repeat id -> (chrom, start, end, motifs, vc_filter or None))
    """
    clusters = {}
    repeats = {}
    with open_maybe_gzipped(path) as f:
        for line in f:
            chrom, start, end, info_field = line.rstrip("\n").split("\t")[:4]
            locus_id = locus_id_of(info_field)
            motifs = motifs_of(info_field)
            if is_variation_cluster(info_field):
                struc = re.search(r"STRUC=(.*)$", info_field).group(1)
                clusters[locus_id] = (chrom, start, end, motifs, struc)
            else:
                repeats[locus_id] = (chrom, start, end, motifs, vc_filter_of(info_field))
    return clusters, repeats


def check_clusters_match(clusters_a, clusters_b, label_a, label_b):
    failures = []
    only_in_a = set(clusters_a) - set(clusters_b)
    only_in_b = set(clusters_b) - set(clusters_a)
    if only_in_a:
        failures.append(f"{len(only_in_a):,d} variation clusters are in {label_a} but not "
                         f"{label_b}, e.g. {sorted(only_in_a)[:3]}")
    if only_in_b:
        failures.append(f"{len(only_in_b):,d} variation clusters are in {label_b} but not "
                         f"{label_a}, e.g. {sorted(only_in_b)[:3]}")

    mismatched = [locus_id for locus_id in clusters_a.keys() & clusters_b.keys()
                  if clusters_a[locus_id] != clusters_b[locus_id]]
    if mismatched:
        example = mismatched[0]
        failures.append(f"{len(mismatched):,d} variation clusters differ between {label_a} and "
                         f"{label_b}, e.g. {example}: {clusters_a[example]} vs {clusters_b[example]}")
    return failures


def check_isolated_repeats_match_expectation(clusters_a, repeats_a, repeats_b, label_a, label_b):
    """Returns failures if the legacy catalog's standalone repeat rows don't match what the
    all_repeats catalog implies they should be.
    """
    member_ids = {member_id for cluster in clusters_a.values() for member_id in
                  members_of_cluster(f"STRUC={cluster[4]}")}

    flagged_ids = {locus_id for locus_id, fields in repeats_a.items() if fields[4]}
    expected_isolated_ids = set(repeats_a) - flagged_ids - member_ids
    actual_isolated_ids = set(repeats_b)

    failures = []
    missing = expected_isolated_ids - actual_isolated_ids
    if missing:
        failures.append(f"{len(missing):,d} repeats are unflagged, non-member rows in {label_a} "
                         f"but have no row in {label_b}, e.g. {sorted(missing)[:3]}")

    extra = actual_isolated_ids - expected_isolated_ids
    # Distinguish the two ways a legacy row can be unexpected, since they point at different bugs.
    extra_members = extra & member_ids
    extra_non_members = extra - member_ids

    if extra_members:
        failures.append(f"{len(extra_members):,d} rows in {label_b} are variation cluster "
                         f"members in {label_a}, so they should not also have a standalone row, "
                         f"e.g. {sorted(extra_members)[:3]}")

    unexpected = extra_non_members - LOCI_ABSENT_FROM_REPEAT_CATALOG
    if unexpected:
        failures.append(f"{len(unexpected):,d} rows in {label_b} have no corresponding "
                         f"unflagged row in {label_a} at all, e.g. {sorted(unexpected)[:3]}")

    stale = LOCI_ABSENT_FROM_REPEAT_CATALOG - extra_non_members
    if stale:
        failures.append(f"{len(stale):,d} ids in LOCI_ABSENT_FROM_REPEAT_CATALOG are no longer "
                         f"unmatched rows in {label_b}, so they are stale and must be removed from "
                         f"that list, e.g. {sorted(stale)[:3]}")

    return failures, len(member_ids), len(flagged_ids)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("variation_clusters_and_all_repeats_bed_path",
                        help="Path of the variation_clusters_and_all_repeats TRGT catalog")
    parser.add_argument("legacy_variation_clusters_and_isolated_trs_bed_path",
                        help="Path of the legacy variation_clusters_and_isolated_TRs TRGT catalog")
    args = parser.parse_args()

    label_a = "the all_repeats catalog"
    label_b = "the legacy isolated_TRs catalog"

    clusters_a, repeats_a = read_catalog(args.variation_clusters_and_all_repeats_bed_path)
    clusters_b, repeats_b = read_catalog(args.legacy_variation_clusters_and_isolated_trs_bed_path)

    isolated_failures, member_count, flagged_count = check_isolated_repeats_match_expectation(
        clusters_a, repeats_a, repeats_b, label_a, label_b)

    print(f"{label_a}: {len(clusters_a):,d} clusters, {len(repeats_a):,d} repeat rows "
          f"({member_count:,d} cluster members, {flagged_count:,d} flagged VC_FILTER)")
    print(f"{label_b}: {len(clusters_b):,d} clusters, {len(repeats_b):,d} repeat rows")
    print(f"Exempting {len(LOCI_ABSENT_FROM_REPEAT_CATALOG):,d} loci that are in the variation "
          f"clusters TSV but not in the repeat catalog, so only {label_b} can hold a row for them")

    failures = check_clusters_match(clusters_a, clusters_b, label_a, label_b) + isolated_failures

    for failure in failures:
        print(f"FAILED: {failure}")
    if failures:
        raise SystemExit(1)
    print("All checks passed")


if __name__ == "__main__":
    main()
