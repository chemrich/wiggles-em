"""Tests for the port abstraction — the whole coupling surface to PyMOL.

The BridgePort tests run against a **real TCP server** speaking MCPymol's
protocol, in a thread. That verifies the wire layer for real: framing,
half-close, chunking, error mapping. It does not verify PyMOL semantics —
nothing here knows what a protein is — but it means BridgePort is not merely
asserted against a mock of itself.
"""

from __future__ import annotations

import json
import socket
import threading

import pytest

from wiggles_em.port import (
    COMMAND_ACTIONS,
    DATA_ACTIONS,
    BridgePort,
    FakePort,
    PortError,
    PymolPort,
    SendRequestPort,
    interpret,
)

# ── FakePort ────────────────────────────────────────────────────────────────


class TestUnknownActionsAreRefused:
    """FakePort used to answer "OK" to any action outside DATA_ACTIONS.

    That is the gap MCPymol's PR #58 found in its own live sweep: a fake that
    says yes to everything cannot tell a real command from a misspelled one,
    nor from one the viewer would reject on the *meaning* of its arguments
    rather than its signature. Every test in this package passes explicit
    arguments through such a fake, so the blind spot covered the whole suite.
    """

    def test_an_allowlisted_command_is_still_fire_and_forget(self):
        """Commands whose result nothing reads need no stub, or every test
        would carry a wall of noise that hides the queries that matter."""
        port = FakePort()

        assert port.query("isosurface", "surf", "map", 1.0) == "OK"
        assert port.query("delete", "obj") == "OK"

    def test_a_misspelled_action_raises_rather_than_answering_ok(self):
        """One character wrong is the case the old behaviour could not see."""
        port = FakePort()

        with pytest.raises(KeyError, match="does not recognise"):
            port.query("isosurfce", "surf", "map", 1.0)

    def test_an_action_nobody_declared_raises(self):
        port = FakePort()

        with pytest.raises(KeyError, match="does not recognise"):
            port.query("reinitialize")

    def test_a_data_action_without_a_stub_still_raises_its_own_error(self):
        """The two failures are distinct: a missing stub is a hole in the
        test, an unknown action is a bug in the caller."""
        port = FakePort()

        with pytest.raises(KeyError, match="no response for data action"):
            port.query("count_states", "obj")

    def test_get_is_a_data_action_because_its_answer_is_parsed(self):
        """`normalisation_state` parses what `get` returns. While `get` was
        absent from DATA_ACTIONS an unstubbed call answered "OK", which parses
        to None — "PyMOL will not say" — so a test that forgot to stub it was
        indistinguishable from a session that genuinely could not answer.
        """
        assert "get" in DATA_ACTIONS

        with pytest.raises(KeyError, match="no response for data action"):
            FakePort().query("get", "normalize_ccp4_maps")

    def test_the_two_action_sets_do_not_overlap(self):
        """An action in both would take the DATA_ACTIONS branch and the
        allowlist entry would be a lie about what the caller does with it."""
        assert not (DATA_ACTIONS & COMMAND_ACTIONS)

    def test_a_stubbed_response_wins_over_both_sets(self):
        """A test may stub a command to observe it, and that must keep
        working — the allowlist decides what needs no stub, not what may
        not have one."""
        port = FakePort({"delete": "gone"})

        assert port.query("delete", "obj") == "gone"


def test_every_action_the_package_issues_is_declared():
    """The allowlist is only a guard while it matches what the code calls.

    Scans the source for `call(port, "...")` and `port.query("...")` and
    requires each name to be declared. A new action added without a decision
    about which set it belongs in fails here rather than silently answering
    "OK" — which is the failure mode this whole mechanism exists to remove.
    """
    import re
    from pathlib import Path

    src = Path(__file__).resolve().parent.parent / "src" / "wiggles_em"
    pattern = re.compile(
        r'(?:call\(\s*(?:self\.)?port,\s*|(?:port|self\.port)\.query\(\s*)"([a-z_0-9]+)"'
    )

    found: dict[str, str] = {}
    for path in src.rglob("*.py"):
        if path.name == "port.py":
            continue  # defines the sets; its own examples are not call sites
        for name in pattern.findall(path.read_text()):
            found.setdefault(name, str(path.relative_to(src)))

    assert found, "the scan found no call sites at all, so it proves nothing"

    undeclared = {n: where for n, where in found.items() if n not in DATA_ACTIONS | COMMAND_ACTIONS}
    assert not undeclared, (
        f"actions issued but not declared in DATA_ACTIONS or COMMAND_ACTIONS: {undeclared}"
    )


def test_both_ports_satisfy_the_protocol():
    """If this breaks, one adapter may not be substitutable for the other."""
    assert isinstance(FakePort(), PymolPort)
    assert isinstance(BridgePort(), PymolPort)


def test_commands_are_recorded_in_order():
    port = FakePort()
    port.do("one")
    port.do("two")
    assert port.commands == ["one", "two"]
    assert port.transcript == "one\ntwo"


def test_queries_are_keyed_by_action():
    port = FakePort({"count_states": 4})
    assert port.query("count_states", "obj") == 4
    assert port.queried("count_states")
    assert not port.queried("get_coords")


def test_query_arguments_are_recorded():
    port = FakePort({"get_coords": []})
    port.query("get_coords", "obj", 3, state=7)
    assert port.queries == [("get_coords", ("obj", 3), {"state": 7})]


def test_callable_response_receives_the_arguments():
    """How a test varies an answer by state without threading state through."""
    port = FakePort({"get_coords": lambda _sel, state: [(float(state), 0.0, 0.0)]})
    assert port.query("get_coords", "obj", 2) == [(2.0, 0.0, 0.0)]


def test_unstubbed_data_action_raises_rather_than_returning_none():
    """A view asking for data the test did not set up must fail loudly —
    returning None would render an empty view silently."""
    port = FakePort({"count_states": 2})
    with pytest.raises(KeyError, match="no response for data action"):
        port.query("iterate_to_list", "obj", "q")


def test_unstubbed_command_returns_ok_and_is_still_recorded():
    """Commands are fire-and-forget: their result is not consumed, so stubbing
    every one would be noise that hides the data queries that matter."""
    port = FakePort()
    assert port.query("show", "sticks", "obj") == "OK"
    assert port.called("show", "sticks", "obj")


def test_unstubbed_data_action_error_lists_known_actions():
    port = FakePort({"alpha": 1, "beta": 2})
    with pytest.raises(KeyError) as exc:
        port.query("get_names", "objects")
    assert "alpha" in str(exc.value) and "beta" in str(exc.value)


def test_responses_are_copied_not_aliased():
    responses = {"count_states": 1}
    port = FakePort(responses)
    responses["get_names"] = ["x"]
    with pytest.raises(KeyError):
        port.query("get_names", "objects")


def test_helpers():
    port = FakePort()
    port.do("show sticks, obj")
    port.do("color red, obj")
    assert port.ran("sticks") and not port.ran("cartoon")
    assert port.commands_matching("obj") == ["show sticks, obj", "color red, obj"]


# ── BridgePort against a real socket ────────────────────────────────────────


class FakePlugin:
    """A localhost server speaking MCPymol's request/response framing."""

    def __init__(self, handler):
        self.handler = handler
        self.requests: list[dict] = []
        self._sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self._sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self._sock.bind(("127.0.0.1", 0))
        self._sock.listen(1)
        self.port = self._sock.getsockname()[1]
        self._thread = threading.Thread(target=self._serve, daemon=True)
        self._thread.start()

    def _serve(self):
        try:
            conn, _ = self._sock.accept()
        except OSError:
            return
        with conn:
            chunks = []
            while True:
                chunk = conn.recv(65536)
                if not chunk:
                    break
                chunks.append(chunk)
                try:
                    json.loads(b"".join(chunks).decode())
                except json.JSONDecodeError:
                    continue
                break
            raw = b"".join(chunks)
            if not raw:
                return
            request = json.loads(raw.decode())
            self.requests.append(request)
            reply = self.handler(request)
            if reply is not None:
                conn.sendall(reply if isinstance(reply, bytes) else json.dumps(reply).encode())
            try:
                conn.shutdown(socket.SHUT_WR)
            except OSError:
                pass

    def close(self):
        self._sock.close()


@pytest.fixture
def plugin():
    made = []

    def _make(handler):
        server = FakePlugin(handler)
        made.append(server)
        return server, BridgePort(host="127.0.0.1", port=server.port, timeout=5.0)

    yield _make
    for server in made:
        server.close()


def test_do_sends_the_command_and_returns_the_result(plugin):
    server, port = plugin(lambda r: {"status": "success", "result": "Executed"})

    assert port.do("show cartoon, obj") == "Executed"
    assert server.requests[0] == {
        "action": "do",
        "args": ["show cartoon, obj"],
        "kwargs": {},
    }


def test_query_sends_action_args_and_kwargs(plugin):
    server, port = plugin(lambda r: {"status": "success", "result": [1, 2]})

    assert port.query("get_coords", "obj", state=2) == [1, 2]
    assert server.requests[0] == {
        "action": "get_coords",
        "args": ["obj"],
        "kwargs": {"state": 2},
    }


def test_iterate_to_list_round_trips_atom_rows(plugin):
    rows = [["A", "1", "MET", "CA", "", 1.0, 20.0]]
    _, port = plugin(lambda r: {"status": "success", "result": rows})

    assert port.query("iterate_to_list", "polymer", "chain, resi") == rows


def test_plugin_error_becomes_a_PortError(plugin):
    _, port = plugin(lambda r: {"status": "error", "error": "Invalid selection"})

    with pytest.raises(PortError, match="Invalid selection"):
        port.query("count_states", "obj")


def test_missing_iterate_to_list_says_where_to_get_it(plugin):
    """The one action stock MCPymol lacks — the error should not be a puzzle."""
    _, port = plugin(
        lambda r: {
            "status": "error",
            "error": "Unknown action or method not found on cmd: iterate_to_list",
        }
    )

    with pytest.raises(PortError, match="per-atom properties"):
        port.query("iterate_to_list", "obj", "q")


def test_connection_refused_is_actionable():
    """Nothing listening: say what to check, not just 'connection refused'."""
    port = BridgePort(host="127.0.0.1", port=1, timeout=2.0)
    with pytest.raises(PortError, match="Is PyMOL running"):
        port.do("show cartoon")


def test_empty_response_is_an_error(plugin):
    _, port = plugin(lambda r: None)  # accept, reply with nothing
    with pytest.raises(PortError, match="empty response"):
        port.do("noop")


def test_malformed_response_is_an_error(plugin):
    _, port = plugin(lambda r: b"{not json")
    with pytest.raises(PortError, match="malformed response"):
        port.do("noop")


def test_non_object_response_is_an_error(plugin):
    _, port = plugin(lambda r: b"[1, 2, 3]")
    with pytest.raises(PortError, match="expected a JSON object"):
        port.do("noop")


def test_response_split_across_chunks_is_reassembled(plugin):
    """Guards the drain loop: a reply arriving in pieces must still parse."""
    big = ["x" * 1000] * 100
    _, port = plugin(lambda r: {"status": "success", "result": big})

    assert port.query("get_names") == big


def test_do_tolerates_a_result_free_success(plugin):
    _, port = plugin(lambda r: {"status": "success"})
    assert port.do("noop") == "OK"


# ── SendRequestPort — the in-MCPymol adapter ────────────────────────────────


def _responder(**by_action):
    """Build a send_request-shaped callable from {action: response}."""
    calls = []

    def send_request(action, args=None, kwargs=None, **extra):
        calls.append((action, list(args or []), dict(kwargs or {})))
        if action in by_action:
            return by_action[action]
        return {"status": "success", "result": "OK"}

    send_request.calls = calls  # type: ignore[attr-defined]
    return send_request


def test_sendrequestport_satisfies_the_protocol():
    assert isinstance(SendRequestPort(_responder()), PymolPort)


def test_sendrequestport_forwards_action_args_and_kwargs():
    send = _responder(get_coords={"status": "success", "result": [[1.0, 2.0, 3.0]]})
    port = SendRequestPort(send)

    assert port.query("get_coords", "obj", state=2) == [[1.0, 2.0, 3.0]]
    assert send.calls[0] == ("get_coords", ["obj"], {"state": 2})


def test_sendrequestport_turns_an_error_dict_into_an_exception():
    """MCPymol returns error dicts; Wiggles' views are written against
    exceptions, and a dict that looks like a result is the failure mode."""
    send = _responder(count_states={"status": "error", "error": "Invalid selection"})
    with pytest.raises(PortError, match="Invalid selection"):
        SendRequestPort(send).query("count_states", "obj")


def test_sendrequestport_reports_a_missing_iterate_to_list_usefully():
    send = _responder(
        iterate_to_list={"status": "error", "error": "Unknown action or method not found"}
    )
    with pytest.raises(PortError, match="per-atom properties"):
        SendRequestPort(send).query("iterate_to_list", "obj", "q")


def test_sendrequestport_wraps_transport_failures():
    def boom(action, args=None, kwargs=None, **extra):
        raise ConnectionRefusedError("nothing listening")

    with pytest.raises(PortError, match="could not be sent"):
        SendRequestPort(boom).do("show cartoon")


def test_both_adapters_interpret_a_response_identically():
    """The seam's point: one translation, shared. If these diverge, behaviour
    depends on how Wiggles happens to be installed."""
    ok = {"status": "success", "result": 42}
    assert interpret("count_states", ok) == SendRequestPort(_responder(count_states=ok)).query(
        "count_states", "obj"
    )
