"""A gzip writer whose compression runs on a thread (zlib releases the GIL) while the caller goes on reading."""
import queue
import threading
import zlib


class GzipWriter:
    """Writes bytes to a gzip file at level 1; a thread compresses while the caller reads and filters."""

    def __init__(self, path):
        self.file = open(path, 'wb')
        self.compressor = zlib.compressobj(1, zlib.DEFLATED, 31)        # 31: a gzip header
        self.chunks = queue.Queue(maxsize=8)
        self.pending = []
        self.pending_size = 0
        self.error = None
        self.thread = threading.Thread(target=self._compress, daemon=True)
        self.thread.start()

    def _compress(self):
        try:
            while True:
                data = self.chunks.get()
                if data is None:
                    break
                self.file.write(self.compressor.compress(data))
            self.file.write(self.compressor.flush())
            self.file.close()
        except BaseException as e:
            self.error = e

    def write(self, data):
        if len(data) > 1 << 22:
            # in pieces, so that the queue holds at most 8 x 4 MB beyond this one object
            self._hand_over()
            view = memoryview(data)
            for i in range(0, len(data), 1 << 22):
                self.chunks.put(view[i:i + (1 << 22)])
            return
        self.pending.append(data)
        self.pending_size += len(data)
        if self.pending_size >= 1 << 22:
            self._hand_over()

    def _hand_over(self):
        if self.pending:
            self.chunks.put(b''.join(self.pending))
            self.pending = []
            self.pending_size = 0

    def close(self):
        self._hand_over()
        self.chunks.put(None)
        self.thread.join()
        if self.error is not None:
            raise self.error
